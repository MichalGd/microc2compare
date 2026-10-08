"""Dependency-free configuration, metadata and provenance contracts."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import os
import re
import tempfile
from pathlib import Path

STAGES = ("preflight", "register", "qc", "pool", "compartments", "insulation",
          "loops", "consensus", "tads", "quantify", "aggregates", "differential",
          "depth_sensitivity", "report")
SWITCHES = {"qc": "RUN_QC", "pool": "RUN_POOL", "compartments": "RUN_COMPARTMENTS",
            "insulation": "RUN_INSULATION", "tads": "RUN_TADS", "loops": "RUN_LOOPS", "aggregates": "RUN_AGGREGATES",
            "differential": "RUN_DIFFERENTIAL", "depth_sensitivity": "RUN_DEPTH_SENSITIVITY",
            "report": "RUN_FINAL_REPORT"}
PATH_KEYS = {"SAMPLESHEET", "TECHNICAL_MEMBERSHIP", "COMPARISONS", "REFERENCE_REGISTRY",
             "INDEPENDENT_FEATURES_REGISTRY", "OUTDIR", "SCRATCH_ROOT", "CORE_ENV_PREFIX",
             "HICREP_ENV_PREFIX", "DIFF_R_ENV_PREFIX", "POOLS_MANIFEST", "SADDLE_REFERENCE_TRACK",
             "CALIBRATION_REGISTRY", "CORE_LOCK", "HICREP_LOCK", "DIFF_R_LOCK"}
IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
SAMPLE_FIELDS = ("sample_id", "project", "condition", "assay", "biological_replicate",
                 "subject_id", "reference_id", "matrix_path", "provenance_path")


class ContractError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise ContractError(message)


def native_path(value):
    if os.name == "nt" and re.match(r"^/[A-Za-z]/", str(value)):
        value = value[1] + ":" + value[2:]
    return Path(value)


def resolve(value, base):
    if not value:
        return ""
    path, sep, group = str(value).partition("::")
    p = native_path(path)
    if not p.is_absolute():
        p = native_path(base) / p
    return str(p.resolve()) + (sep + group if sep else "")


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def file_hash(path):
    result = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def fingerprint(path, previous=None):
    p = native_path(str(path).split("::")[0])
    require(p.is_file(), f"Missing readable input: {p}")
    stat = p.stat()
    identity = {"path": str(p.resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if previous and all(previous.get(k) == v for k, v in identity.items()):
        identity["sha256"] = previous["sha256"]
    else:
        identity["sha256"] = file_hash(p)
    return identity


def atomic_json(path, data):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    fd, temp = tempfile.mkstemp(prefix=".pending-", dir=p.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(data, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp, p)
    finally:
        if os.path.exists(temp):
            os.unlink(temp)


def read_json(path):
    with open(path, encoding="utf-8-sig") as handle:
        return json.load(handle)


def read_table(path, fields=()):
    p = Path(path)
    require(p.is_file(), f"Missing manifest: {p}")
    with p.open(encoding="utf-8-sig", newline="") as handle:
        first = handle.readline()
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter="\t" if "\t" in first else ",")
        require(reader.fieldnames and len(reader.fieldnames) == len(set(reader.fieldnames)),
                f"Missing/duplicate headers in {p}")
        require(set(fields) <= set(reader.fieldnames), f"{p}: required headers {list(fields)}")
        rows = list(reader)
    for line, row in enumerate(rows, 2):
        require(None not in row and all(v is not None for v in row.values()),
                f"{p}:{line}: row has the wrong number of fields")
    return rows


def write_table(path, rows, fields):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_structured(path):
    if str(path).endswith(".json"):
        return read_json(path)
    import yaml
    with open(path, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def normalize(environment, config_dir):
    keys = read_json(Path(__file__).parents[2] / "schemas/config.schema.json")["properties"]
    c = {}
    for key, definition in keys.items():
        value = environment.get(key, definition.get("default"))
        if definition["type"] == "boolean":
            require(str(value).lower() in ("true", "false"), f"{key}: use true or false")
            value = str(value).lower() == "true"
        elif definition["type"] == "integer":
            require(str(value).isdigit(), f"{key}: expected nonnegative integer")
            value = int(value)
            require(value >= definition.get("minimum", 0), f"{key}: below minimum")
            require(value <= definition.get("maximum", math.inf), f"{key}: above maximum")
        elif definition["type"] == "number":
            try:
                value = float(value)
            except (TypeError, ValueError) as error:
                raise ContractError(f"{key}: expected a number") from error
            require(math.isfinite(value) and definition.get("minimum", -math.inf) <= value <=
                    definition.get("maximum", math.inf), f"{key}: outside valid range")
        else:
            value = "" if value is None else str(value)
        require("enum" not in definition or value in definition["enum"],
                f"{key}: permitted values {definition.get('enum')}")
        if key in PATH_KEYS:
            value = resolve(value, config_dir)
        c[key] = value
    for key in ("PROJECT", "ANALYSIS_SET_ID", "RUN_ID", "REFERENCE_ID"):
        require(IDENTIFIER.fullmatch(c[key]), f"{key}: supply a stable safe identifier")
    require(c["SCHEMA_VERSION"] == 1, "Only schema version 1 is supported")
    require(c["PIPELINE_NAME"]=="microc2compare","PIPELINE_NAME must remain microc2compare")
    require(c["LOOP_COUNT_WINDOW_BINS"]==1,"v0.1 supports one raw anchor-pair pixel, not wider count windows")
    require(c["FORMAL_MIN_UNITS_PER_CONDITION"]>=3 and c["PRELIMINARY_MIN_UNITS_PER_CONDITION"]>=2 and c["MIN_RESIDUAL_DF"]>=2,
            "Statistical policy floors cannot be weakened below 3 formal/2 preliminary units and 2 residual df")
    c["run_dir"] = str(Path(c["OUTDIR"]) / c["PROJECT"] / c["ANALYSIS_SET_ID"] / c["RUN_ID"])
    c["config_dir"] = str(native_path(config_dir).resolve())
    from .analyses import configure_selection
    configure_selection(c, environment.get("M2C_ANALYSES", ""))
    require(not c["RUN_POOL"] or c["POOL_MODE"] != "import" or c["POOLS_MANIFEST"],
            "RUN_POOL=true with POOL_MODE=import requires POOLS_MANIFEST")
    require(not (c["RUN_DEPTH_SENSITIVITY"] and c["DEPTH_SENSITIVITY_MODE"] == "balanced_features")
            or c["BALANCE_NEW_DERIVED_MAPS"], "Balanced depth sensitivity requires BALANCE_NEW_DERIVED_MAPS=true")
    require(not c["RUN_AGGREGATES"] or "saddle" not in c["AGGREGATE_MODULES"].split(",") or c["SADDLE_PARTITION_SOURCE"] != "all_sample_pool" or c["RUN_POOL"],
            "Pool-based saddle partition requires RUN_POOL=true and a validated pool")
    require(not c["RUN_AGGREGATES"] or "saddle" not in c["AGGREGATE_MODULES"].split(",") or c["SADDLE_PARTITION_SOURCE"] != "external" or c["SADDLE_REFERENCE_TRACK"],
            "External saddle partition requires SADDLE_REFERENCE_TRACK")
    require(not c["LOOP_SENSITIVITY_RESOLUTIONS"] or c["LOOP_SENSITIVITY_RESOLUTIONS"] == "5000",
            "Only optional 5000-bp loop sensitivity is supported")
    require(not c["LOOP_SENSITIVITY_RESOLUTIONS"] or c["FINE_RESOLUTION_QC_APPROVED"],
            "5-kb calling requires FINE_RESOLUTION_QC_APPROVED=true after pilot review")
    windows = [int(x) for x in c["INSULATION_WINDOWS_BP"].split(",")]
    require(c["INSULATION_PRIMARY_WINDOW_BP"] in windows and all(x > 0 and x % c["INSULATION_RESOLUTION"] == 0 for x in windows),
            "Insulation windows must include the primary window and be positive bin multiples")
    modules = c["AGGREGATE_MODULES"].split(",")
    require(set(modules) <= {"apa", "ada", "saddle"} and len(set(modules)) == len(modules), "AGGREGATE_MODULES: use distinct apa,ada,saddle names")
    require(c["TAD_MAX_LENGTH_BP"] >= c["TAD_MIN_LENGTH_BP"], "TAD maximum length is below minimum")
    thresholds = [float(x) for x in c["TAD_SENSITIVITY_THRESHOLDS"].split(",") if x]
    require(all(math.isfinite(x) and x >= 0 for x in thresholds), "TAD sensitivity thresholds must be finite and nonnegative")
    if c["INSULATION_BOUNDARY_THRESHOLD"] not in ("Li", "Otsu"):
        try: threshold = float(c["INSULATION_BOUNDARY_THRESHOLD"])
        except ValueError as error: raise ContractError("Boundary threshold must be Li, Otsu or a nonnegative number") from error
        require(math.isfinite(threshold) and threshold >= 0, "Boundary threshold must be finite and nonnegative")
    return c


def enabled(c, stage):
    from .analyses import needs_loop_counts
    if stage=="qc":return c["RUN_QC"] or c["RUN_HICREP"]
    if stage == "consensus":
        return c["RUN_LOOPS"] or c["RUN_INSULATION"]
    if stage == "quantify":
        return needs_loop_counts(c)
    if stage == "tads":
        return c.get("RUN_TADS", False) and c["RUN_INSULATION"]
    return c.get(SWITCHES.get(stage), True)


def allocation_limits():
    """Effective affinity/cgroup allocation; configured budgets are never enlarged."""
    cpu = float(len(os.sched_getaffinity(0))) if hasattr(os, "sched_getaffinity") else float(os.cpu_count() or 1)
    memory = None
    if platform.system() == "Linux":
        info = Path("/proc/meminfo").read_text()
        match = re.search(r"^MemAvailable:\s+(\d+)", info, re.M)
        if match: memory = int(match[1]) * 1024
        # Inspect every ancestor: a parent cgroup can impose a tighter limit.
        entries = Path("/proc/self/cgroup").read_text().splitlines()
        group = next((x.split(":",2)[2] for x in entries if x.startswith("0::")), None)
        if group is not None:
            base = Path("/sys/fs/cgroup")
            current = base / group.lstrip("/")
            while current.is_relative_to(base):
                quota = current / "cpu.max"
                if quota.exists():
                    amount, period = quota.read_text().split()
                    if amount != "max": cpu = min(cpu, int(amount)/int(period))
                limit = current / "memory.max"
                if limit.exists() and limit.read_text().strip() != "max":
                    remaining = int(limit.read_text())
                    usage = current / "memory.current"
                    if usage.exists(): remaining -= int(usage.read_text())
                    memory = remaining if memory is None else min(memory, remaining)
                if current == base: break
                current = current.parent
        elif any(x.split(":")[1] for x in entries):
            raise ContractError("Legacy cgroup v1 allocation is not supported; use cgroup v2 or a reviewed unrestricted host")
    return {"effective_cpus":cpu,"available_memory_bytes":memory}


def validate_allocation(c):
    limits = allocation_limits()
    require(c["TOTAL_CPUS"] <= limits["effective_cpus"],
            f"TOTAL_CPUS={c['TOTAL_CPUS']} exceeds effective affinity/cgroup allocation {limits['effective_cpus']}; lower the budget")
    if limits["available_memory_bytes"] is not None:
        require(c["MAX_MEMORY_MB"]*1000000 <= limits["available_memory_bytes"],
                "MAX_MEMORY_MB exceeds available physical/cgroup memory; preserve headroom and lower the budget")
    return limits


def validate_metadata(c):
    samples = read_table(c["SAMPLESHEET"], SAMPLE_FIELDS)
    require(samples, "The biological sample manifest is empty")
    ids = [s["sample_id"] for s in samples]
    require(len(set(ids)) == len(ids), "sample_id must be unique")
    base = Path(c["SAMPLESHEET"]).parent
    subjects = set()
    for s in samples:
        for key in SAMPLE_FIELDS:
            require(s[key], f"{s.get('sample_id')}: empty {key}")
        require(IDENTIFIER.fullmatch(s["sample_id"]), f"Unsafe sample_id: {s['sample_id']}")
        require(s["sample_id"] not in ("mask", "labels") and not s["sample_id"].startswith("pool__"), "sample_id uses a reserved internal name")
        require(s["assay"] in ("microc", "hic"), f"{s['sample_id']}: assay must be microc or hic")
        require(s["reference_id"] == c["REFERENCE_ID"], f"{s['sample_id']}: reference mismatch")
        require(s["project"] == c["PROJECT"], f"{s['sample_id']}: project mismatch")
        pair = (s["subject_id"], s["condition"])
        require(pair not in subjects, f"Repeated subject/condition {pair}: v0.1 requires one consolidated observation")
        subjects.add(pair)
        for field in ("matrix_path", "provenance_path", "pairs_path", "pair_stats_path"):
            if s.get(field):
                s[field] = resolve(s[field], base)
    matrices = [s["matrix_path"] for s in samples]
    require(len(set(matrices)) == len(matrices), "Same matrix cannot supply two biological observations")
    require(c["ANALYSIS_MODE"] == "assay_comparison" or len({s["assay"] for s in samples}) == 1,
            "Mixed assays require explicit descriptive ANALYSIS_MODE=assay_comparison")
    membership = read_table(c["TECHNICAL_MEMBERSHIP"],
                            ("sample_id", "library_id", "source_id", "source_kind", "dedup_group_id", "dedup_scope", "source_checksum"))
    owners, sources, libraries = {}, set(), {}
    for row in membership:
        require(row["sample_id"] in ids, f"Unknown membership sample: {row['sample_id']}")
        require(all(row.values()), f"Incomplete membership: {row}")
        require(row["source_id"] not in sources, f"Duplicate source/constituent: {row['source_id']}")
        sources.add(row["source_id"])
        library = row["library_id"]
        require(library not in owners or owners[library] == row["sample_id"], f"Library {library} belongs to multiple observations")
        owners[library] = row["sample_id"]
        libraries.setdefault(library, []).append(row)
        require(row["source_kind"] in ("lane", "resequencing", "independent_library"), "Unknown source_kind")
    require(set(owners.values()) == set(ids), "Every biological sample needs technical-library membership")
    for library, rows in libraries.items():
        if len(rows) > 1 and any(r["source_kind"] in ("lane", "resequencing") for r in rows):
            require(all(r["dedup_scope"] == "joint_library" for r in rows) and
                    len({r["dedup_group_id"] for r in rows}) == 1,
                    f"{library}: split-library rows require upstream joint library deduplication")
    comparisons = read_table(c["COMPARISONS"], ("comparison_id", "case_condition", "reference_condition", "design", "inference", "modules")) if c["RUN_DIFFERENTIAL"] or enabled(c, "tads") else []
    if c.get("_comparison_modules") and not enabled(c, "tads"):
        comparisons = [dict(row, modules=";".join(m for m in row["modules"].split(";") if m in c["_comparison_modules"])) for row in comparisons]
        comparisons = [row for row in comparisons if row["modules"]]
    require(len({x["comparison_id"] for x in comparisons}) == len(comparisons), "Duplicate comparison IDs")
    for row in comparisons:
        require(IDENTIFIER.fullmatch(row["comparison_id"]), "Unsafe comparison ID")
        require(row["inference"] in ("descriptive", "preliminary", "formal"), "Unknown inference mode")
        require(row["case_condition"] != row["reference_condition"], "Contrast conditions must differ")
        active_comparison = c["RUN_DIFFERENTIAL"] and bool(set(row["modules"].split(";")) & set(c.get("_comparison_modules", ["loops", "insulation"])))
        terms = formula_terms(row["design"]) if active_comparison else ["condition"]
        require("condition" in terms, "Design must include condition")
        chosen = [s for s in samples if s["condition"] in (row["case_condition"], row["reference_condition"])]
        require({s["condition"] for s in chosen} == {row["case_condition"], row["reference_condition"]}, "Comparison missing a condition")
        for s in chosen:
            require(all(s.get(term) for term in terms), f"{s['sample_id']}: missing design covariate")
        require(set(row["modules"].split(";")) <= {"loops", "insulation"}, "Formal compartments/domains/stripes/pixels are deferred")
        if active_comparison and row["inference"] != "descriptive":
            require(c["ANALYSIS_MODE"] != "assay_comparison", "Cross-assay inference is deferred")
            minimum = c["FORMAL_MIN_UNITS_PER_CONDITION"] if row["inference"] == "formal" else c["PRELIMINARY_MIN_UNITS_PER_CONDITION"]
            for condition in (row["case_condition"], row["reference_condition"]):
                require(len({s["subject_id"] for s in chosen if s["condition"] == condition}) >= minimum,
                        f"{row['comparison_id']}: {condition} requires {minimum} independent units")
            for s in chosen:
                require(all(r["dedup_scope"] != "unknown" for r in membership if r["sample_id"] == s["sample_id"]), "Unknown dedup scope blocks inference")
            if len({s["subject_id"] for s in chosen}) < len(chosen):
                require("subject_id" in terms, "Repeated biological subjects require subject_id in design")
    return {"samples": samples, "comparisons": comparisons, "membership": membership}


def formula_terms(formula):
    require(re.fullmatch(r"\s*~\s*[A-Za-z_]\w*(?:\s*\+\s*[A-Za-z_]\w*)*\s*", formula),
            "Only additive ~column+column formulas are supported; no executable expressions")
    terms = [x.strip() for x in formula.strip().lstrip("~").split("+")]
    require(len(set(terms)) == len(terms), "Duplicate design terms")
    return terms


def concurrency(c, threads, memory, matrix=False, io=False):
    require(threads > 0 and memory > 0, "Task reservations must be positive")
    count = min(c["MAX_PARALLEL_JOBS"], c["MAX_PARALLEL_SAMPLES"],
                c["TOTAL_CPUS"] // threads, c["MAX_MEMORY_MB"] // memory)
    if matrix:
        count = min(count, c["MAX_PARALLEL_MATRIX"])
    if io:
        count = min(count, c["MAX_PARALLEL_IO"])
    require(count >= 1, f"Task needs {threads} CPUs/{memory} MB, exceeding the configured budget")
    return count
