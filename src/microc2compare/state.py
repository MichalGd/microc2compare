"""Fixed-stage signatures, output validation and resource admission."""
from __future__ import annotations

import importlib.metadata
from contextlib import contextmanager
from contextvars import ContextVar
from functools import wraps
import json
import os
import platform
import time
from pathlib import Path

from .contracts import (STAGES, ContractError, atomic_json, concurrency, digest, enabled,
                        file_hash, fingerprint, load_structured, resolve, read_json, require, validate_metadata)

RESERVATIONS = {"preflight": "PREFLIGHT", "register": "PREFLIGHT", "qc": "QC",
                "pool": "POOL", "compartments": "COMPARTMENTS", "insulation": "INSULATION",
                "loops": "LOOPS", "consensus": "QUANTIFY", "tads": "TADS", "quantify": "QUANTIFY",
                "chromosight": "CHROMOSIGHT", "chromosight_discovery": "CHROMOSIGHT", "annotation": "QUANTIFY",
                "stripes": "STRIPES", "stripe_consensus": "PREFLIGHT", "browser_export": "BROWSER_EXPORT",
                "aggregates": "AGGREGATES", "differential": "DIFFERENTIAL",
                "depth_sensitivity": "DEPTH_SENSITIVITY", "report": "REPORT", "hicrep": "HICREP"}
SAMPLE_STAGES = {"qc", "compartments", "insulation", "loops", "chromosight", "chromosight_discovery", "stripes"}


class _ValidationAudit:
    """One synchronous dependency audit; never retained between public calls."""

    def __init__(self, config):
        self.config = config
        self.active = True
        self.results = {"complete": {}, "signature": {}, "scopes": {}}
        self.source_digests = {}
        self.source_hashes = {}
        self.metadata = None


_VALIDATION_AUDIT = ContextVar("microc2compare_validation_audit", default=None)


def _current_validation_audit():
    audit = _VALIDATION_AUDIT.get()
    # Copied execution contexts may still reference an audit that has finished.
    return audit if audit is not None and audit.active else None


@contextmanager
def _validation_context(c):
    audit = _current_validation_audit()
    if audit is not None and audit.config is c:
        yield audit
        return
    # A nested audit of another config must not borrow identities from this run.
    audit = _ValidationAudit(c)
    token = _VALIDATION_AUDIT.set(audit)
    try:
        yield audit
    finally:
        audit.active = False
        _VALIDATION_AUDIT.reset(token)


def _cached_validation(kind):
    def decorate(function):
        @wraps(function)
        def wrapped(c, *args, **kwargs):
            with _validation_context(c) as audit:
                # Keeping the uncached body usable also makes exact equivalence
                # measurable in regression tests without a runtime cache switch.
                if audit is None:
                    return function(c, *args, **kwargs)
                cache = audit.results[kind]
                key = (args, tuple(sorted(kwargs.items())))
                if key not in cache:
                    cache[key] = function(c, *args, **kwargs)
                result = cache[key]
                return list(result) if kind == "scopes" else result
        return wrapped
    return decorate


def _source_file_hash(path):
    audit = _current_validation_audit()
    if audit is None:
        return file_hash(path)
    key = Path(path)
    if key not in audit.source_hashes:
        audit.source_hashes[key] = file_hash(path)
    return audit.source_hashes[key]


def source_digest(stage):
    audit = _current_validation_audit()
    if audit is not None and stage in audit.source_digests:
        return audit.source_digests[stage]
    root = Path(__file__).parents[2]
    names = ["contracts.py", "state.py", "analyses.py", "matrix.py", "chromosomes.py", "resolutions.py", "features.py", "fine_loops.py", "fine_validation.py", "tooling.py", "annotation.py", "science.py"]
    if stage == "tads": names.append("tads.py")
    if stage in ("stripes", "stripe_consensus"): names.append("stripes.py")
    if stage in ("loops", "chromosight_discovery", "chromosight"): names.append("loop_profiles.py")
    if stage == "aggregates": names += ["aggregation.py", "expected_cache.py", "parallel.py"]
    if stage == "quantify": names += ["quantification.py", "parallel.py"]
    if stage == "qc": names += ["qc_parallel.py", "expected_cache.py", "parallel.py"]
    if stage == "browser_export": names += [p.name for p in (root / "src/microc2compare").glob("browser_*.py")]
    if stage == "report":
        names = ["contracts.py", "state.py", "analyses.py", "report.py", "fine_report.py", "stripe_report.py", "volcano.py", "resolutions.py", "loop_profiles.py"]
    sources = {n: _source_file_hash(root / "src/microc2compare" / n) for n in names if (root / "src/microc2compare" / n).is_file()}
    if stage == "differential":
        sources.update({p.name: _source_file_hash(p) for p in (root / "scripts").glob("differential*.R")})
        sources["calibrate.py"]=_source_file_hash(root/"scripts/calibrate.py")
        sources["calibration-protocol.md"]=_source_file_hash(root/"validation/calibration-protocol.md")
    sources["scientific.py"] = _source_file_hash(root / "scripts/scientific.py")
    for p in [*(root/"schemas").glob("*.json"),*(root/"scripts/lib").glob("*.sh"),root/"scripts/run_stage.sh"]:
        sources[str(p.relative_to(root))]=_source_file_hash(p)
    result = digest(sources)
    if audit is not None:
        audit.source_digests[stage] = result
    return result


def inventory():
    packages = {}
    for name in ("numpy", "scipy", "pandas", "cooler", "cooltools", "bioframe", "hicrep", "h5py", "matplotlib", "PyYAML", "jsonschema", "mustache-hic", "chromosight", "stripepy-hic", "pyBigWig"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    return {"python": platform.python_version(), "platform": platform.platform(), "packages": packages}


@_cached_validation("scopes")
def scopes(c, stage):
    if stage not in SAMPLE_STAGES:
        return ["run"]
    audit = _current_validation_audit()
    if audit is None:
        metadata = validate_metadata(c)
    else:
        if audit.metadata is None:
            audit.metadata = validate_metadata(c)
        metadata = audit.metadata
    ids = [s["sample_id"] for s in metadata["samples"]]
    if stage in ("compartments", "insulation", "loops") and c["RUN_POOL"]:
        registry = Path(c["run_dir"]) / "02_pools/registry.json"
        require(registry.is_file(), "Pool registry unavailable; run enabled pool stage first")
        ids += ["pool__" + p["pool_id"] for p in read_json(registry)["pools"]]
    return ids


def prerequisites(c, stage):
    from .analyses import aggregate_modules, comparison_modules, has_loop_discovery
    if stage == "preflight":
        return []
    if stage == "register":
        return ["preflight"]
    deps = ["register"]
    if stage in ("loops", "chromosight_discovery"):
        deps += ["qc"]
    if stage in ("compartments", "insulation", "loops") and c["RUN_POOL"]:
        deps += ["pool"]
    if stage == "consensus":
        deps += [s for s in ("loops", "chromosight_discovery", "insulation") if enabled(c, s)]
    if stage == "tads": deps += ["insulation"]
    if stage == "quantify" and has_loop_discovery(c):
        deps += ["consensus"]
    if stage in ("chromosight", "annotation"): deps += ["quantify"]
    if stage == "aggregates":
        modules = aggregate_modules(c)
        if "apa" in modules and has_loop_discovery(c): deps += ["consensus"]
        if "ada" in modules: deps += ["insulation", "consensus"]
        if "saddle" in modules: deps += ["compartments"]
    if stage == "differential":
        modules = comparison_modules(c)
        if "loops" in modules: deps += ["quantify"]
        if "insulation" in modules: deps += ["insulation"]
    if stage == "depth_sensitivity":
        deps += ["qc"]
        if c["DEPTH_SENSITIVITY_MODE"] == "balanced_features":
            deps += ["quantify"] + [s for s in ("insulation", "compartments") if enabled(c, s)]
    if stage == "stripe_consensus": deps += ["stripes"]
    if stage == "browser_export":
        deps += [s for s in STAGES[2:STAGES.index("browser_export")] if enabled(c, s)]
    if stage == "report":
        deps += [s for s in STAGES[2:-1] if enabled(c, s)]
        if c["RUN_HICREP"]:deps.append("hicrep")
    return list(dict.fromkeys(deps))


def record_path(c, stage, scope):
    return Path(c["run_dir"]) / "completion" / scope / (stage + ".done.json")


def relevant_config(c, stage):
    common = ("REFERENCE_ID", "ANALYSIS_MODE", "WEIGHT_COLUMN", "WEIGHT_CONVENTION", "LEGACY_WEIGHT_VALIDATION", "CHROMOSOME_SCOPE", "SEX_CHROMOSOMES", "ANALYSIS_CHROMOSOMES", "MATRIX_REGISTRY")
    prefixes = {
        "preflight": tuple(c),
        "register": ("COMPARTMENT_RES", "INSULATION_RES", "LOOP_RES", "LOOP_SENSITIVITY", "PS_RES", "HICREP_RES"),
        "qc": ("PS_", "INSULATION_IGNORE", "LOOP_MIN_DISTANCE", "LOOP_MAX_DISTANCE"),
        "hicrep": ("HICREP_",), "pool": ("POOL_", "BALANCE_NEW", "KEEP_POOLED"),
        "compartments": ("COMPARTMENT_", "MIN_PHASING", "SADDLE_"),
        "insulation": ("INSULATION_",), "loops": ("LOOP_", "MUSTACHE_", "FINE_RESOLUTION", "RESOLUTION_ELIGIBILITY", "HSIEH_METHODS"),
        "tads": ("TAD_", "INSULATION_"),
        "consensus": ("LOOP_CONSENSUS", "LOOP_DISCOVERY", "LOOP_RESOLUTION", "CROSS_RESOLUTION", "INSULATION_", "RUN_LOOPS", "RUN_INSULATION"),
        "chromosight": ("CHROMOSIGHT_", "LOOP_"), "annotation": ("ANNOTATION_",),
        "chromosight_discovery": ("CHROMOSIGHT_", "LOOP_", "RESOLUTION_ELIGIBILITY", "MAX_DENSE_BINS"),
        "stripes": ("STRIPE_",), "stripe_consensus": ("STRIPE_",),
        "browser_export": ("BROWSER_", "EXPORT_", "HIGLASS_"),
        "quantify": ("LOOP_", "INDEPENDENT_FEATURES"), "aggregates": ("SADDLE_", "LOOP_", "INSULATION_", "AGGREGATE_MODULES"),
        "differential": ("LOOP_", "INSULATION_", "FDR_", "FORMAL_", "PRELIMINARY_", "MIN_RESIDUAL_", "CALIBRATION_", "ALLOW_UNCALIBRATED"),
        "depth_sensitivity": ("DEPTH_", "DOWNSAMPLING_", "RANDOM_SEED", "BALANCE_NEW", "KEEP_DOWNSAMPLED", "LOOP_", "INSULATION_"),
        "report": tuple(c),
    }[stage]
    result = {k: v for k, v in c.items() if k in common or any(k.startswith(p) for p in prefixes)}
    if stage == "register":
        from .analyses import needs_loop_counts, resolution_set
        result = {k: c[k] for k in common}
        result["resolutions"] = resolution_set(c)
        result["balanced_needed"] = any(c[k] for k in ("RUN_QC", "RUN_COMPARTMENTS", "RUN_INSULATION", "RUN_LOOPS", "RUN_AGGREGATES", "RUN_STRIPES")) or needs_loop_counts(c)
    if stage == "differential": result["comparison_modules"] = c.get("_comparison_modules", ["loops", "insulation"])
    if stage == "aggregates":
        from .analyses import aggregate_modules
        modules = aggregate_modules(c)
        selected = ["AGGREGATE_MODULES", "INSULATION_IGNORE_DIAGS"]
        if "saddle" in modules: selected += [k for k in c if k.startswith("SADDLE_")]
        if "ada" in modules: selected += ["INSULATION_PRIMARY_WINDOW_BP", "ADA_RESOLUTION"]
        if "apa" in modules: selected += ["LOOP_MIN_DISTANCE_BP", "LOOP_MAX_DISTANCE_BP", "INDEPENDENT_FEATURES_REGISTRY", "LOOP_RESOLUTION", "LOOP_QUANTIFICATION_RESOLUTIONS", "APA_FLANK_BP"]
        result = {k: c[k] for k in set(common) | set(selected)}
    return result


@_cached_validation("signature")
def signature(c, stage, scope):
    root = Path(c["run_dir"])
    input_file = root / "00_metadata/inputs.json"
    inputs = read_json(input_file) if input_file.exists() else {}
    # Verify immutable inputs against the registered identity, without repeated full hashing.
    for entry in inputs.values():
        p = Path(entry["path"])
        require(p.is_file() and p.stat().st_size == entry["size"] and p.stat().st_mtime_ns == entry["mtime_ns"],
                f"Input changed since registration: {p}; restart registration")
    deps = {}
    for dependency in prerequisites(c, stage):
        require(enabled(c, dependency), f"{stage}: required {dependency} stage disabled")
        for ds in scopes(c, dependency):
            p = record_path(c, dependency, ds)
            require(p.exists(), f"{stage}: missing current {dependency}/{ds}; run prerequisite")
            record = read_json(p)
            require(complete(c, dependency, ds), f"{stage}: stale/invalid {dependency}/{ds}; resume prerequisite")
            # Preflight admission always refreshes. Parameters unrelated to bin
            # registration must not invalidate all scientific measurements.
            if dependency != "preflight": deps[f"{dependency}/{ds}"] = record["signature"]
    envfile = root / "00_metadata/environments.json"
    environments = read_json(envfile) if envfile.exists() else {}
    if stage not in ("preflight", "report", "differential"):
        required={"hicrep"} if stage=="hicrep" else {"core", "browser"} if stage == "browser_export" else {"core", "stripes"} if stage in ("stripes", "stripe_consensus") else {"core", "loops"} if stage in ("loops","chromosight", "chromosight_discovery") else {"core"}
        environments = {k: v for k, v in environments.items() if k in required}
    # Keep immutable identities checked above, but bind optional assets only to
    # consumers. A comparison/threshold experiment must not redo unrelated QC.
    irrelevant = set()
    for key, consumers in (("COMPARISONS", {"preflight", "tads", "differential", "report"}),
                           ("CALIBRATION_REGISTRY", {"preflight", "differential", "report"}),
                           ("INDEPENDENT_FEATURES_REGISTRY", {"preflight", "quantify", "aggregates", "differential", "report"}),
                           ("SADDLE_REFERENCE_TRACK", {"preflight", "aggregates", "report"}),
                           ("POOLS_MANIFEST", {"preflight", "pool", "report"})):
        if stage not in consumers and c.get(key): irrelevant.add(c[key])
    inputs = {k: v for k, v in inputs.items() if k not in irrelevant}
    return digest({"stage": stage, "scope": scope, "inputs": inputs,
                   "config": relevant_config(c, stage), "source": source_digest(stage),
                   "environments": environments, "dependencies": deps})


@_cached_validation("complete")
def complete(c, stage, scope):
    p = record_path(c, stage, scope)
    if not p.exists():
        return False
    try:
        r = read_json(p)
        if r.get("status") not in ("success", "success_empty") or r["signature"] != signature(c, stage, scope):
            return False
        for output in r["outputs"]:
            q = Path(c["run_dir"]) / output["path"]
            if output.get("discarded_intermediate"):
                continue
            if not q.is_file() or q.stat().st_size != output["size"] or file_hash(q) != output["sha256"]:
                return False
        return True
    except (OSError, KeyError, ValueError, RecursionError):
        return False


def publish_completion(c, stage, scope, outputs, empty=False, extra=None):
    root = Path(c["run_dir"])
    entries = []
    for p in outputs:
        p = Path(p)
        require(p.is_file(), f"Missing required output {p}")
        if p.resolve().is_relative_to(root.resolve()):
            output_path = str(p.resolve().relative_to(root.resolve()))
        else:
            require(stage == "browser_export" and c.get("BROWSER_OUTPUT_DIR") and p.resolve().is_relative_to(Path(c["BROWSER_OUTPUT_DIR"]).resolve()), "Output escapes the run-owned directory")
            output_path = str(p.resolve())
        entries.append({"path": output_path, "size": p.stat().st_size, "sha256": file_hash(p)})
    atomic_json(record_path(c, stage, scope), {"stage": stage, "scope": scope,
        "status": "success_empty" if empty else "success", "signature": signature(c, stage, scope),
        "outputs": entries, "completed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "provenance": extra or {}})


def stage_concurrency(c, stage):
    key = RESERVATIONS[stage]
    threads, memory = c["THREADS_" + key], c["MEMORY_MB_" + key]
    if stage in ("pool", "depth_sensitivity") and c["BALANCE_NEW_DERIVED_MAPS"]:
        threads=max(threads,c["THREADS_BALANCE"])
        memory=max(memory,c["MEMORY_MB_BALANCE"])
    count = concurrency(c, threads, memory, matrix=stage in SAMPLE_STAGES or stage in ("pool", "aggregates", "depth_sensitivity", "hicrep"),
                        io=stage in ("register", "pool", "depth_sensitivity"))
    if c.get("LIMIT_MATRIX_READERS", False) and stage in SAMPLE_STAGES:
        # Conservatively count one potential reader per reserved worker. Unlike
        # a job cap, this also bounds concurrent readers inside observations.
        require(c["MAX_PARALLEL_IO"] >= threads, "Matrix-reader cap is below this stage's worker reservation; lower its threads or raise MAX_PARALLEL_IO")
        count = min(count, c["MAX_PARALLEL_IO"] // threads)
    return count


def restore_required_pools(c, stage):
    """Rebuild discarded explicit new pools only if an unfinished consumer needs them."""
    if stage not in ("compartments","insulation","loops","aggregates") or not c["RUN_POOL"]:
        return
    registry = Path(c["run_dir"])/"02_pools/registry.json"
    if not registry.exists() or all(complete(c,stage,s) for s in scopes(c,stage)):
        return
    pools = read_json(registry)["pools"]
    missing = any(p.get("task_owned") and any(not Path(u.split("::")[0]).exists() for u in p["uris"].values()) for p in pools)
    if missing:
        require(c["POOL_MODE"] == "create" and c["BALANCE_NEW_DERIVED_MAPS"],
                "Discarded derived pools need explicit create/balance settings before their consumers can restart")
        record_path(c,"pool","run").unlink(missing_ok=True)
        return True
    return False


def scope_identity(c):
    """Freeze effective chromosome lengths AND selected view intervals without Cooler."""
    from .chromosomes import chromosome_selection
    registry=load_structured(c['REFERENCE_REGISTRY']);base=Path(c['REFERENCE_REGISTRY']).parent
    sizes=[line.split() for line in Path(resolve(registry['chromsizes'],base)).read_text().splitlines() if line.strip() and not line.startswith('#')]
    views=[line.split() for line in Path(resolve(registry['views'],base)).read_text().splitlines() if line.strip() and not line.startswith('#')]
    selection=chromosome_selection(c,registry,[r[0] for r in sizes],[r[0] for r in views])
    selected=set(selection['chromosomes'])
    return dict(selection,coordinate_digest=digest({'sizes':[r for r in sizes if r[0] in selected], 'views':[r for r in views if r[0] in selected]}))


def initialize(c, pid, hostname):
    root = Path(c["run_dir"])
    # Reject scope migration BEFORE changing any saved metadata or completion product.
    config_path=root/"00_metadata/config.json"
    policy_keys=("CHROMOSOME_SCOPE","SEX_CHROMOSOMES","ANALYSIS_CHROMOSOMES")
    if config_path.exists():
        previous=read_json(config_path)
        if all(key in previous for key in policy_keys):
            require(all(previous[key]==c[key] for key in policy_keys),
                    "Chromosome scope changed: choose a new RUN_ID; previous results must be preserved")
        else:
            require(False,"Legacy run has no frozen chromosome policy: use its original release to resume or choose a new RUN_ID")
    identity=scope_identity(c)
    scope_path=root/'00_metadata/chromosome_policy.json'
    if scope_path.exists():
        require(read_json(scope_path)==identity,'Effective chromosome coordinates/scope changed: use a new RUN_ID')
    metadata = validate_metadata(c)
    old_path = root / "00_metadata/inputs.json"
    old = read_json(old_path) if old_path.exists() else {}
    keys = ["SAMPLESHEET", "TECHNICAL_MEMBERSHIP", "REFERENCE_REGISTRY"]
    if c["RUN_DIFFERENTIAL"] or enabled(c, "tads"): keys.append("COMPARISONS")
    from .analyses import aggregate_modules, needs_loop_counts
    if needs_loop_counts(c) or "apa" in aggregate_modules(c): keys.append("INDEPENDENT_FEATURES_REGISTRY")
    if c["RUN_POOL"]: keys.append("POOLS_MANIFEST")
    if "saddle" in aggregate_modules(c): keys.append("SADDLE_REFERENCE_TRACK")
    for key in ("MATRIX_REGISTRY","RESOLUTION_ELIGIBILITY_REGISTRY","HSIEH_METHODS_REGISTRY","ANNOTATION_REGISTRY"):
        if c.get(key): keys.append(key)
    paths = [c[k] for k in keys if c[k]]
    paths += [s[k] for s in metadata["samples"] for k in ("matrix_path", "provenance_path")]
    if c.get("MATRIX_REGISTRY"):
        from .resolutions import sample_matrix
        from .analyses import resolution_set
        paths += [sample_matrix(c,s,res) for s in metadata["samples"] for res in resolution_set(c)]
    if c["RUN_DIFFERENTIAL"] and any(x["inference"]!="descriptive" for x in metadata["comparisons"]) and Path(c["CALIBRATION_REGISTRY"]).is_file():
        paths.append(c["CALIBRATION_REGISTRY"])
    # Keep previously registered reference/feature/pool assets in restart identity checks.
    paths += [p for p in old if p not in paths and Path(p).is_file()]
    inputs = {str(p).split("::")[0]: fingerprint(p, old.get(str(p).split("::")[0])) for p in paths}
    atomic_json(old_path, inputs)
    atomic_json(scope_path,identity)
    atomic_json(root / "00_metadata/config.json", c)
    atomic_json(root / "00_metadata/samples.json", metadata)
    atomic_json(root / ".run-lock/owner.json", {"pid": pid, "hostname": hostname, "run_id": c["RUN_ID"]})


def invalidate(c, stage):
    affected = {stage}
    if stage=="qc" and c["RUN_HICREP"]:affected.add("hicrep")
    for following in STAGES:
        if set(prerequisites(c, following)) & affected:
            affected.add(following)
    for p in (Path(c["run_dir"]) / "completion").glob("*/*.done.json"):
        if p.name.removesuffix(".done.json") in affected:
            p.unlink()
    return sorted(affected)
