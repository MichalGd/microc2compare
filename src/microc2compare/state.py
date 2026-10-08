"""Fixed-stage signatures, output validation and resource admission."""
from __future__ import annotations

import importlib.metadata
import json
import os
import platform
import time
from pathlib import Path

from .contracts import (STAGES, ContractError, atomic_json, concurrency, digest, enabled,
                        file_hash, fingerprint, read_json, require, validate_metadata)

RESERVATIONS = {"preflight": "PREFLIGHT", "register": "PREFLIGHT", "qc": "QC",
                "pool": "POOL", "compartments": "COMPARTMENTS", "insulation": "INSULATION",
                "loops": "LOOPS", "consensus": "QUANTIFY", "tads": "QUANTIFY", "quantify": "QUANTIFY",
                "aggregates": "AGGREGATES", "differential": "DIFFERENTIAL",
                "depth_sensitivity": "DEPTH_SENSITIVITY", "report": "REPORT", "hicrep": "HICREP"}
SAMPLE_STAGES = {"qc", "compartments", "insulation", "loops"}


def source_digest(stage):
    root = Path(__file__).parents[2]
    names = ["contracts.py", "state.py", "analyses.py", "matrix.py", "features.py", "science.py"]
    if stage == "tads": names.append("tads.py")
    if stage == "report":
        names = ["contracts.py", "state.py", "analyses.py", "report.py", "volcano.py"]
    sources = {n: file_hash(root / "src/microc2compare" / n) for n in names if (root / "src/microc2compare" / n).is_file()}
    if stage == "differential":
        sources.update({p.name: file_hash(p) for p in (root / "scripts").glob("differential*.R")})
        sources["calibrate.py"]=file_hash(root/"scripts/calibrate.py")
        sources["calibration-protocol.md"]=file_hash(root/"validation/calibration-protocol.md")
    sources["scientific.py"] = file_hash(root / "scripts/scientific.py")
    for p in [*(root/"schemas").glob("*.json"),*(root/"scripts/lib").glob("*.sh"),root/"scripts/run_stage.sh"]:
        sources[str(p.relative_to(root))]=file_hash(p)
    return digest(sources)


def inventory():
    packages = {}
    for name in ("numpy", "scipy", "pandas", "cooler", "cooltools", "bioframe", "hicrep", "h5py", "matplotlib", "PyYAML", "jsonschema"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    return {"python": platform.python_version(), "platform": platform.platform(), "packages": packages}


def scopes(c, stage):
    if stage not in SAMPLE_STAGES:
        return ["run"]
    ids = [s["sample_id"] for s in validate_metadata(c)["samples"]]
    if stage in ("compartments", "insulation", "loops") and c["RUN_POOL"]:
        registry = Path(c["run_dir"]) / "02_pools/registry.json"
        require(registry.is_file(), "Pool registry unavailable; run enabled pool stage first")
        ids += ["pool__" + p["pool_id"] for p in read_json(registry)["pools"]]
    return ids


def prerequisites(c, stage):
    from .analyses import aggregate_modules, comparison_modules
    if stage == "preflight":
        return []
    if stage == "register":
        return ["preflight"]
    deps = ["register"]
    if stage == "loops":
        deps += ["qc"]
    if stage in ("compartments", "insulation", "loops") and c["RUN_POOL"]:
        deps += ["pool"]
    if stage == "consensus":
        deps += [s for s in ("loops", "insulation") if enabled(c, s)]
    if stage == "tads": deps += ["insulation"]
    if stage == "quantify" and c["RUN_LOOPS"]:
        deps += ["consensus"]
    if stage == "aggregates":
        modules = aggregate_modules(c)
        if "apa" in modules and c["RUN_LOOPS"]: deps += ["consensus"]
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
    if stage == "report":
        deps += [s for s in STAGES[2:-1] if enabled(c, s)]
        if c["RUN_HICREP"]:deps.append("hicrep")
    return list(dict.fromkeys(deps))


def record_path(c, stage, scope):
    return Path(c["run_dir"]) / "completion" / scope / (stage + ".done.json")


def relevant_config(c, stage):
    common = ("REFERENCE_ID", "ANALYSIS_MODE", "WEIGHT_COLUMN", "WEIGHT_CONVENTION", "LEGACY_WEIGHT_VALIDATION")
    prefixes = {
        "preflight": tuple(c),
        "register": ("COMPARTMENT_RES", "INSULATION_RES", "LOOP_RES", "LOOP_SENSITIVITY", "PS_RES", "HICREP_RES"),
        "qc": ("PS_", "INSULATION_IGNORE", "LOOP_MIN_DISTANCE", "LOOP_MAX_DISTANCE"),
        "hicrep": ("HICREP_",), "pool": ("POOL_", "BALANCE_NEW", "KEEP_POOLED"),
        "compartments": ("COMPARTMENT_", "MIN_PHASING", "SADDLE_"),
        "insulation": ("INSULATION_",), "loops": ("LOOP_", "FINE_RESOLUTION"),
        "tads": ("TAD_", "INSULATION_"),
        "consensus": ("LOOP_CONSENSUS", "INSULATION_", "RUN_LOOPS", "RUN_INSULATION"),
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
        result["balanced_needed"] = any(c[k] for k in ("RUN_QC", "RUN_COMPARTMENTS", "RUN_INSULATION", "RUN_LOOPS", "RUN_AGGREGATES")) or needs_loop_counts(c)
    if stage == "differential": result["comparison_modules"] = c.get("_comparison_modules", ["loops", "insulation"])
    if stage == "aggregates":
        from .analyses import aggregate_modules
        modules = aggregate_modules(c)
        selected = ["AGGREGATE_MODULES", "INSULATION_IGNORE_DIAGS", "LOOP_RESOLUTION"]
        if "saddle" in modules: selected += [k for k in c if k.startswith("SADDLE_")]
        if "ada" in modules: selected += ["INSULATION_PRIMARY_WINDOW_BP"]
        if "apa" in modules: selected += ["LOOP_MIN_DISTANCE_BP", "LOOP_MAX_DISTANCE_BP", "INDEPENDENT_FEATURES_REGISTRY"]
        result = {k: c[k] for k in set(common) | set(selected)}
    return result


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
        environments = {k: v for k, v in environments.items() if k == ("hicrep" if stage == "hicrep" else "core")}
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
        entries.append({"path": str(p.relative_to(root)), "size": p.stat().st_size, "sha256": file_hash(p)})
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
    return concurrency(c, threads, memory, matrix=stage in SAMPLE_STAGES or stage in ("pool", "aggregates", "depth_sensitivity", "hicrep"),
                       io=stage in ("register", "pool", "depth_sensitivity"))


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


def initialize(c, pid, hostname):
    root = Path(c["run_dir"])
    metadata = validate_metadata(c)
    old_path = root / "00_metadata/inputs.json"
    old = read_json(old_path) if old_path.exists() else {}
    keys = ["SAMPLESHEET", "TECHNICAL_MEMBERSHIP", "REFERENCE_REGISTRY"]
    if c["RUN_DIFFERENTIAL"] or enabled(c, "tads"): keys.append("COMPARISONS")
    from .analyses import aggregate_modules, needs_loop_counts
    if needs_loop_counts(c) or "apa" in aggregate_modules(c): keys.append("INDEPENDENT_FEATURES_REGISTRY")
    if c["RUN_POOL"]: keys.append("POOLS_MANIFEST")
    if "saddle" in aggregate_modules(c): keys.append("SADDLE_REFERENCE_TRACK")
    paths = [c[k] for k in keys if c[k]]
    paths += [s[k] for s in metadata["samples"] for k in ("matrix_path", "provenance_path")]
    if c["RUN_DIFFERENTIAL"] and any(x["inference"]!="descriptive" for x in metadata["comparisons"]) and Path(c["CALIBRATION_REGISTRY"]).is_file():
        paths.append(c["CALIBRATION_REGISTRY"])
    # Keep previously registered reference/feature/pool assets in restart identity checks.
    paths += [p for p in old if p not in paths and Path(p).is_file()]
    inputs = {str(p).split("::")[0]: fingerprint(p, old.get(str(p).split("::")[0])) for p in paths}
    atomic_json(old_path, inputs)
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
