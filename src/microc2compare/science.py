"""Scientific stage implementations using existing maps and explicit derived tasks."""
from __future__ import annotations

import itertools
import heapq
import json
import os
import shutil
import subprocess
from pathlib import Path

import cooler
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from .contracts import (IDENTIFIER, atomic_json, digest, file_hash, fingerprint, formula_terms,
                        load_structured, read_json, read_table, require, resolve,
                        validate_metadata, validate_allocation, write_table)
from .features import (LOOP_COLUMNS, boundary_consensus, canonical_loops,
                       consensus_loops, feature_id, quantify, quantify_many)
from .matrix import (MaskedCooler, aligned_pixel_chunks, analysis_bins, bounded_raw, domain, expected_table, matrix_uri,
                     native_domain, native_indices, pixel_chunks, raw_summary, reference,
                     resolution_set, scope_record, tool_views, unordered_total, view_labels, weights)
from .state import inventory, source_digest, scopes as stage_scopes
from .analyses import aggregate_modules, comparison_modules, has_loop_discovery, needs_loop_counts
from .resolutions import loop_resolutions, quantification_resolutions, sample_matrix


def save(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary=path.with_name("."+path.name+".pending")
    frame.to_csv(temporary, sep="\t", index=False, na_rep="NA")
    os.replace(temporary,path)
    return path


def track(frame,column,path):
    finite=frame.loc[np.isfinite(frame[column]),["chrom","start","end",column]]
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    pending=path.with_name("."+path.name+".pending")
    finite.to_csv(pending,sep="\t",index=False,header=False)
    os.replace(pending,path)
    # An entirely unavailable track is represented by a comment, never by zeros.
    if path.stat().st_size==0:path.write_text("# No finite eligible measurements\n")
    return path


def design_matrix(samples, comparison, c):
    selected = [s for s in samples if s["condition"] in (comparison["reference_condition"], comparison["case_condition"])]
    columns, names, contrast = [np.ones(len(selected))], ["intercept"], [0.0]
    for term in formula_terms(comparison["design"]):
        values = [s[term] for s in selected]
        if term not in ("condition","subject_id","batch","assay"):
            try: numeric = np.array(values, dtype=float)
            except ValueError: numeric = None
        else: numeric = None
        if numeric is not None:
            require(np.isfinite(numeric).all(), f"Nonfinite covariate {term}")
            columns.append(numeric); names.append(term); contrast.append(0.0)
        else:
            levels = sorted(set(values))
            if term == "condition": levels = [comparison["reference_condition"], comparison["case_condition"]]
            for level in levels[1:]:
                columns.append(np.array([float(v == level) for v in values]))
                names.append(term + "=" + level)
                contrast.append(1.0 if term == "condition" and level == comparison["case_condition"] else 0.0)
    x = np.column_stack(columns)
    if comparison["inference"] != "descriptive":
        require(np.linalg.matrix_rank(x) == x.shape[1], f"{comparison['comparison_id']}: rank-deficient/confounded design")
        require(len(selected)-x.shape[1] >= c["MIN_RESIDUAL_DF"], f"{comparison['comparison_id']}: insufficient residual degrees of freedom")
    return {"sample_ids": [s["sample_id"] for s in selected], "columns": names,
            "matrix": x.tolist(), "contrast": contrast, "residual_df": int(len(selected)-np.linalg.matrix_rank(x))}


def external_features(c, samples):
    if not c["INDEPENDENT_FEATURES_REGISTRY"]:
        return pd.DataFrame(columns=LOOP_COLUMNS), None
    registry = load_structured(c["INDEPENDENT_FEATURES_REGISTRY"])
    require(isinstance(registry, dict), "Feature registry must describe one frozen BEDPE set")
    import jsonschema
    jsonschema.validate(registry, read_json(Path(__file__).parents[2]/"schemas/features.schema.json"))
    require(registry["reference_id"] == c["REFERENCE_ID"], "External feature reference mismatch")
    path = resolve(registry["path"], Path(c["INDEPENDENT_FEATURES_REGISTRY"]).parent)
    require(file_hash(path) == registry["sha256"], "Frozen external feature checksum mismatch")
    require(not set(registry["discovery_samples"]) & {s["sample_id"] for s in samples}, "External loop discovery reuses analysis samples")
    require(registry["independent_of_analysis_samples"] is True, "Formal loop features must be independently discovered")
    require(registry.get("source") and registry.get("independence_evidence"), "External candidates require inspectable source/independence evidence")
    features=canonical_loops(pd.read_csv(path, sep="\t"), c["REFERENCE_ID"], registry["resolution"])
    features['representative_resolution']=registry['resolution']
    if Path(c["REFERENCE_REGISTRY"]).is_file():
        selected=reference(c)[0]["effective_scope"]["chromosomes"]
        included=features.chrom1.isin(selected)&features.chrom2.isin(selected)
        registry=dict(registry,out_of_scope_feature_count=int((~included).sum()))
        features=features.loc[included].reset_index(drop=True)
    return features, registry


def environment_record(c, kind):
    if c["ENVIRONMENT_MODE"] == "development":
        if kind in ("core","loops","stripes","browser"): return inventory()
        if kind == "hicrep": return {"hicrep": inventory()["packages"]["hicrep"]}
        command = [c["R_SCRIPT"], "--vanilla", "-e", 'cat(R.version.string,"\\n"); for(p in c("edgeR","limma","statmod")) cat(p,as.character(packageVersion(p)),"\\n")']
    else:
        prefix = c[{"core": "CORE_ENV_PREFIX", "loops":"LOOP_ENV_PREFIX", "hicrep": "HICREP_ENV_PREFIX", "diff_r": "DIFF_R_ENV_PREFIX", "stripes": "STRIPE_ENV_PREFIX", "browser": "BROWSER_ENV_PREFIX"}[kind]]
        if kind=="loops" and not prefix: prefix=c["CORE_ENV_PREFIX"]
        require((Path(prefix)/"conda-meta").is_dir(), f"Missing {kind} prefix {prefix}; run explicit bootstrap")
        packages = subprocess.run([c["CONDA_EXE"], "list", "-p", prefix, "--json"], check=True, text=True, capture_output=True)
        if c["ENVIRONMENT_MODE"] == "locked":
            lock = c[{"core":"CORE_LOCK","loops":"LOOP_LOCK","hicrep":"HICREP_LOCK","diff_r":"DIFF_R_LOCK","stripes":"STRIPE_LOCK","browser":"BROWSER_LOCK"}[kind]]
            require(lock and Path(lock+".inventory.json").is_file(), f"Locked {kind} mode needs a tested lock inventory")
            expected = read_json(lock+".inventory.json")
            require(digest(json.loads(packages.stdout)) == expected["package_digest"], f"{kind} package drift from lock inventory")
        if kind != "diff_r": return {"packages": json.loads(packages.stdout), "core_runtime": inventory() if kind == "core" else None}
        command = [c["CONDA_EXE"], "run", "--no-capture-output", "-p", prefix, "Rscript", "--vanilla", "-e", 'cat(R.version.string,"\\n"); for(p in c("edgeR","limma","statmod")) cat(p,as.character(packageVersion(p)),"\\n")']
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    return {"runtime": result.stdout.strip(), "packages": json.loads(packages.stdout) if c["ENVIRONMENT_MODE"] != "development" else None}


def calibration_design_domain(c,comparison,samples):
    selected=[s for s in samples if s["condition"] in (comparison["reference_condition"],comparison["case_condition"])]
    design=design_matrix(selected,comparison,c)
    return {"reference_units":len({s["subject_id"] for s in selected if s["condition"]==comparison["reference_condition"]}),
            "case_units":len({s["subject_id"] for s in selected if s["condition"]==comparison["case_condition"]}),
            "observations":len(selected),"subjects":len({s["subject_id"] for s in selected}),
            "columns":len(design["columns"]),"residual_df":int(design["residual_df"])}


def calibration_key(c, module, comparison, environments, design_domain=None):
    keys = ("LOOP_RESOLUTION", "LOOP_MIN_DISTANCE_BP", "LOOP_MAX_DISTANCE_BP", "LOOP_MIN_TOTAL_COUNT", "LOOP_MIN_MEAN_COUNT_AT_MIN_EXPOSURE", "LOOP_MIN_FEATURES_FORMAL", "INSULATION_RESOLUTION", "INSULATION_PRIMARY_WINDOW_BP", "INSULATION_MIN_VALID_PIXEL_FRACTION", "FDR_THRESHOLD", "CHROMOSOME_SCOPE", "SEX_CHROMOSOMES", "ANALYSIS_CHROMOSOMES", "LOOP_COUNT_WINDOW_BINS")
    if design_domain is None:design_domain=calibration_design_domain(c,comparison,validate_metadata(c)["samples"])
    return digest({"source": source_digest("differential"), "environment": environments,"design_domain":design_domain,
                   "module": module, "design": formula_terms(comparison["design"]),
                   "recipe": {k:c[k] for k in keys},"mask_policy":"intersection_of_usable_bins","reference_id":c['REFERENCE_ID']})


def check_calibration(c, comparison, environments):
    if comparison["inference"] == "descriptive": return "descriptive"
    if c["ALLOW_UNCALIBRATED_DEVELOPMENT"]: return "uncalibrated"
    require(Path(c["CALIBRATION_REGISTRY"]).is_file(), "Inference requires a passing calibration registry; run the calibration protocol first")
    artifacts = read_json(c["CALIBRATION_REGISTRY"])["artifacts"]
    for module in comparison["modules"].split(";"):
        key = calibration_key(c, module, comparison, environments)
        match = [x for x in artifacts if x.get("calibration_key") == key and x.get("status") == "passed"]
        require(match, f"{comparison['comparison_id']}/{module}: no matching passing calibration artifact; key={key}")
        for artifact in match:
            evidence = resolve(artifact["results_path"], Path(c["CALIBRATION_REGISTRY"]).parent)
            require(file_hash(evidence) == artifact["results_sha256"], "Calibration results fingerprint mismatch")
            results = read_json(evidence)
            require(results["source_digest"] == source_digest("differential") and results["calibration_key"] == key,
                    "Calibration evidence recipe/source mismatch")
            expected={"equal_depth_null","unequal_depth_null","mixture_10pct","broad_30pct","sparse_null"}
            if "subject_id" in formula_terms(comparison["design"]):expected.add("paired_null")
            scenarios=results.get("primary_scenarios",[])
            require(results.get("status")=="passed" and results.get("module")==module and
                    {s["scenario"] for s in scenarios}==expected and
                    all(s["n_datasets"] >= 200 and s["fdp_upper_95"] <= 0.075 for s in scenarios),
                    "Calibration does not meet prespecified primary acceptance criteria")
    return "calibrated"


def inspect_inputs(c, writing=False):
    import jsonschema
    metadata = validate_metadata(c)
    allocation = validate_allocation(c)
    for target in (c["OUTDIR"],c["SCRATCH_ROOT"]):
        existing=Path(target)
        while not existing.exists():existing=existing.parent
        require(shutil.disk_usage(existing).free>=c["MIN_FREE_DISK_MB"]*1000000,
                f"Insufficient free disk near {target}; require at least MIN_FREE_DISK_MB")
    root = Path(__file__).parents[2]
    for key, schema in (("samples", "samples"), ("membership", "membership"), ("comparisons", "comparisons")):
        jsonschema.validate(metadata[key], read_json(root/"schemas"/(schema+".schema.json")))
    ref, sizes, views = reference(c)
    jsonschema.validate(ref, read_json(root/"schemas/reference.schema.json"))
    if c["RUN_COMPARTMENTS"]:
        phasing_track(c, analysis_bins(sizes, views, c["COMPARTMENT_RESOLUTION"]))
    for s in metadata["samples"]:
        prov = read_json(s["provenance_path"])
        jsonschema.validate(prov, read_json(root/"schemas/provenance.schema.json"))
        own = [r for r in metadata["membership"] if r["sample_id"] == s["sample_id"]]
        require(set(prov["library_ids"]) == {r["library_id"] for r in own} and set(prov["source_ids"]) == {r["source_id"] for r in own}, f"{s['sample_id']}: upstream membership does not match biological manifest")
        if has_loop_discovery(c) or needs_loop_counts(c) or "apa" in aggregate_modules(c):
            require(c["LOOP_MIN_DISTANCE_BP"] >= prov["filters"]["min_distance_bp"], "Loop distance starts below upstream contact cutoff")
    for res in resolution_set(c):
        domain(c, metadata["samples"], res)
    from .fine_validation import validate as validate_fine
    validate_fine(c,metadata['samples'],ref)
    if c["RUN_POOL"] and c["POOLS_MANIFEST"]:
        manifest=read_table(c["POOLS_MANIFEST"])
        jsonschema.validate(manifest,read_json(root/"schemas/pools.schema.json"))
        require(len({p["pool_id"] for p in manifest})==len(manifest),"Duplicate pool identifiers")
        for pool in manifest:
            require(IDENTIFIER.fullmatch(pool["pool_id"]),"Unsafe pool_id; use letters/digits/underscore/dot/hyphen")
            require(pool["pool_scope"] in ("condition","all_samples") and pool["reference_id"]==c["REFERENCE_ID"],"Pool scope/reference mismatch")
            members=pool["member_sample_ids"].split(";")
            wanted={s["sample_id"] for s in metadata["samples"] if pool["pool_scope"]=="all_samples" or s["condition"]==pool["condition"]}
            require(len(members)==len(set(members)) and set(members)==wanted and members,"Pool membership must exactly match its scope")
            for key in ("matrix_path","provenance_path"):
                path=resolve(pool[key],Path(c["POOLS_MANIFEST"]).parent).split("::")[0]
                require(Path(path).is_file(),f"Missing imported pool {key}: {path}")
    feature_set, feature_registry = external_features(c, metadata["samples"]) if "loops" in comparison_modules(c) or "apa" in aggregate_modules(c) or has_loop_discovery(c) else (pd.DataFrame(columns=LOOP_COLUMNS), None)
    environments = {"core": environment_record(c, "core")}
    needed_tools = c["RUN_COMPARTMENTS"] or c["RUN_INSULATION"] or (c["RUN_LOOPS"] and c["LOOP_CALLER"]=="cooltools_dots")
    if needed_tools:
        import cooltools
        require(hasattr(cooltools, "eigs_cis") and hasattr(cooltools, "insulation") and hasattr(cooltools, "dots"), "Core cooltools API incompatibility")
    if c["RUN_HICREP"]:
        environments["hicrep"] = environment_record(c, "hicrep")
    if c["RUN_CHROMOSIGHT"] or c.get("RUN_CHROMOSIGHT_DISCOVERY") or c["RUN_LOOPS"] and c["LOOP_CALLER"]=="mustache":
        environments["loops"]=environment_record(c,"loops")
        from .tooling import probe_command
        tested=subprocess.run(probe_command(c),check=True,capture_output=True,text=True)
        environments["loops"]["adapter_probe"]=json.loads(tested.stdout)
    if c.get("RUN_BROWSER_EXPORT") and c.get("EXPORT_HIGLASS") and c.get("BROWSER_ENV_PREFIX"):
        environments["browser"] = environment_record(c, "browser")
        command = [c["CONDA_EXE"], "run", "--no-capture-output", "-p", c["BROWSER_ENV_PREFIX"], "python"] if c["ENVIRONMENT_MODE"] != "development" else [__import__("sys").executable]
        code = "import sys,json; sys.path.insert(0," + repr(str(root/"src")) + "); from microc2compare.browser_runtime import probe; print(json.dumps(probe()))"
        result = subprocess.run(command + ["-c", code], check=True, capture_output=True, text=True)
        environments["browser"]["probe"] = json.loads(result.stdout)
    if c.get("RUN_STRIPES"):
        environments["stripes"] = environment_record(c, "stripes")
        command = [c["CONDA_EXE"], "run", "--no-capture-output", "-p", c["STRIPE_ENV_PREFIX"], "python"] if c["ENVIRONMENT_MODE"] != "development" else [__import__("sys").executable]
        result = subprocess.run(command + ["-c", "import json; from microc2compare.stripes import probe; print(json.dumps(probe()))"], capture_output=True, text=True)
        # The installed CLI probe is performed in its own prefix; source is supplied explicitly.
        if result.returncode:
            code = "import sys,json; sys.path.insert(0," + repr(str(root/"src")) + "); from microc2compare.stripes import probe; print(json.dumps(probe()))"
            result = subprocess.run(command + ["-c", code], check=True, capture_output=True, text=True)
        environments["stripes"]["probe"] = json.loads(result.stdout)
    if c["RUN_ANNOTATION"]:
        from .annotation import resources
        annotation_resources=resources(c)
    allowed_comparisons = set(c.get("_comparison_modules", ["loops", "insulation"]))
    inferential = c["RUN_DIFFERENTIAL"] and any(x["inference"] != "descriptive" and set(x["modules"].split(";")) & allowed_comparisons for x in metadata["comparisons"])
    if inferential:
        environments["diff_r"] = environment_record(c, "diff_r")
    for comparison in metadata["comparisons"]:
        if c["RUN_DIFFERENTIAL"]:
            if c.get("_comparison_modules"):
                comparison = dict(comparison, modules=";".join(m for m in comparison["modules"].split(";") if m in c["_comparison_modules"]))
                if not comparison["modules"]: continue
            design_matrix(metadata["samples"], comparison, c)
            modules = comparison["modules"].split(";")
            require("insulation" not in modules or c["RUN_INSULATION"], "Insulation comparison requires RUN_INSULATION=true")
            if comparison["inference"] != "descriptive" and "loops" in modules:
                require(feature_registry is not None and len(feature_set) > 0, "Formal loops require independent frozen features")
            if "loops" in modules and comparison["inference"]!="descriptive":
                for res in quantification_resolutions(c):check_calibration(dict(c,LOOP_RESOLUTION=res),comparison,environments)
            else:check_calibration(c, comparison, environments)
    if writing:
        out = Path(c["run_dir"])/"00_metadata"
        prior = read_json(out/"inputs.json")
        for key in ("chromsizes", "views", "blacklist", "phasing_track"):
            if key == "phasing_track" and not c["RUN_COMPARTMENTS"]: continue
            if ref.get(key): prior[ref[key]] = fingerprint(ref[key], prior.get(ref[key]))
        if feature_registry:
            path = resolve(feature_registry["path"], Path(c["INDEPENDENT_FEATURES_REGISTRY"]).parent)
            prior[path] = fingerprint(path, prior.get(path))
        if c["RUN_ANNOTATION"]:
            for key in ("gtf","ccre"):
                if annotation_resources.get(key):
                    path=annotation_resources[key];prior[path]=fingerprint(path,prior.get(path))
        if inferential and Path(c["CALIBRATION_REGISTRY"]).is_file():
            registry_path=c["CALIBRATION_REGISTRY"]
            prior[registry_path]=fingerprint(registry_path,prior.get(registry_path))
            for artifact in read_json(registry_path).get("artifacts",[]):
                if artifact.get("results_path"):
                    evidence=resolve(artifact["results_path"],Path(registry_path).parent)
                    if Path(evidence).is_file():prior[evidence]=fingerprint(evidence,prior.get(evidence))
        if c["RUN_POOL"] and c["POOLS_MANIFEST"]:
            for pool in read_table(c["POOLS_MANIFEST"]):
                for key in ("matrix_path", "provenance_path"):
                    path = resolve(pool[key], Path(c["POOLS_MANIFEST"]).parent)
                    prior[path.split("::")[0]] = fingerprint(path, prior.get(path.split("::")[0]))
        atomic_json(out/"inputs.json", prior)
        atomic_json(out/"environments.json", environments)
        output = out/"preflight.json"
        atomic_json(output, {"status":"validated", "allocation":allocation, "samples":len(metadata["samples"]), "resolutions":resolution_set(c), "validation_status":"uncalibrated" if c["ALLOW_UNCALIBRATED_DEVELOPMENT"] else "software_preflight_only", "normalization":"inherited_upstream_weights"})
        return [output, out/"environments.json"]
    return {"status":"validated", "samples":len(metadata["samples"]), "resolutions":resolution_set(c), "environments":environments}


def register(c, scope):
    samples = validate_metadata(c)["samples"]
    outputs = []
    for res in resolution_set(c):
        bins, common, labels, views, values = domain(c, samples, res)
        path = Path(c["run_dir"])/f"00_metadata/masks/{res}.npz"
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, mask=common, labels=labels, **values)
        outputs += [path, save(bins, path.with_suffix(".bins.tsv")), save(views, path.with_suffix(".views.tsv"))]
        audit = path.with_suffix(".scope.json")
        atomic_json(audit, scope_record(c, samples, res, bins, common))
        outputs.append(audit)
    return outputs


def context(c, scope, resolution):
    base = Path(c["run_dir"])/f"00_metadata/masks/{resolution}"
    bins = pd.read_csv(str(base)+".bins.tsv", sep="\t")
    views = pd.read_csv(str(base)+".views.tsv", sep="\t")
    with np.load(str(base)+".npz") as archive:
        if scope.startswith("pool__"):
            pool = next(p for p in read_json(Path(c["run_dir"])/"02_pools/registry.json")["pools"] if p["pool_id"] == scope[6:])
            uri = pool["uris"][str(resolution)]
            raw = cooler.Cooler(uri)
            pool_config=dict(c,WEIGHT_COLUMN="weight",WEIGHT_CONVENTION="multiplicative") if pool.get("task_owned") else c
            values = weights(raw, pool_config, read_json(pool["provenance_path"]), resolution)[native_indices(raw, bins)]
        else:
            sample = next(s for s in validate_metadata(c)["samples"] if s["sample_id"] == scope)
            uri = matrix_uri(sample_matrix(c, sample, resolution), resolution)
            raw = cooler.Cooler(uri)
            values = archive[scope]
        common = archive["mask"] & np.isfinite(values)
        native, mask, labels, values = native_domain(raw, bins, common, archive["labels"], values)
    return MaskedCooler(uri, values), native, mask, labels, views, values


def qc(c, scope):
    from .qc_parallel import run
    return run(c, scope)


def phasing_track(c, bins):
    ref, _, _ = reference(c)
    require(ref.get("phasing_track"), "Compartments require a registered fixed phasing track")
    track = pd.read_csv(ref["phasing_track"], sep="\t", header=None, names=["chrom", "start", "end", "value"])
    # Exact grid is required; no hidden interpolation/lift-over.
    merged = bins.merge(track, on=["chrom", "start", "end"], how="left", validate="one_to_one")
    merged["value"] = pd.to_numeric(merged["value"], errors="raise")
    require(not np.isinf(merged.value.to_numpy(dtype=float)).any(), "Phasing track contains infinite values")
    require(merged.value.notna().any(), "Phasing track has no matching bins")
    return merged


def registered_table(c, resolution, frame):
    """Export feature measurements on the common ordered grid, not native IDs."""
    base = pd.read_csv(Path(c["run_dir"])/f"00_metadata/masks/{resolution}.bins.tsv", sep="\t")
    frame = frame.copy()
    frame["chrom"] = frame.chrom.astype(str)
    return base.merge(frame, on=["chrom", "start", "end"], how="left", validate="one_to_one")


def compartments(c, scope):
    import cooltools
    clr, bins, mask, labels, views, values = context(c, scope, c["COMPARTMENT_RESOLUTION"])
    phase_values = phasing_track(c, bins)
    dense_sizes=[]
    for v in views.itertuples():
        lo, hi = clr.extent((v.chrom, v.start, v.end))
        require(hi-lo <= c["MAX_DENSE_BINS"], "Compartment view exceeds MAX_DENSE_BINS; increase measured memory/reservation explicitly")
        dense_sizes.append(hi-lo)
    workers=min(c["THREADS_COMPARTMENTS"],len(views))
    estimated=sum(n*n*8*12 for n in sorted(dense_sizes,reverse=True)[:workers])
    require(estimated<=c["MEMORY_MB_COMPARTMENTS"]*1000000,"Concurrent compartment arrays exceed provisional memory bound; reduce workers/views or increase a measured reservation")
    from multiprocess import Pool
    with Pool(workers) as pool:
        eigenvalues, vectors = cooltools.eigs_cis(clr, view_df=tool_views(clr,views), n_eigs=3, clr_weight_name="weight", ignore_diags=c["INSULATION_IGNORE_DIAGS"],map=pool.map)
    vectors = registered_table(c, c["COMPARTMENT_RESOLUTION"], vectors)
    name_column = "name" if "name" in eigenvalues else "region"
    eigenvalues = eigenvalues.set_index(name_column,drop=False).loc[views.name].reset_index(drop=True)
    vectors = vectors.merge(phase_values, on=["chrom", "start", "end"], how="left", validate="one_to_one")
    vectors["oriented_score"] = np.nan
    vectors["selected_component"] = pd.NA
    phasing_qc = []
    for v in views.itertuples():
        selected = (vectors.chrom.astype(str) == v.chrom) & (vectors.start >= v.start) & (vectors.end <= v.end)
        correlations = {}
        for component in ("E1", "E2", "E3"):
            good = selected & np.isfinite(vectors[component]) & np.isfinite(vectors.value)
            if good.sum() >= 3 and vectors.loc[good, component].std() > 0 and vectors.loc[good, "value"].std() > 0:
                correlations[component] = float(pearsonr(vectors.loc[good, component], vectors.loc[good, "value"])[0])
        if correlations:
            winner = max(correlations, key=lambda k: (abs(correlations[k]), -int(k[1:])))
            correlation = correlations[winner]
            if abs(correlation) >= c["MIN_PHASING_CORRELATION"]:
                vectors.loc[selected, "oriented_score"] = vectors.loc[selected, winner] * (1 if correlation >= 0 else -1)
                vectors.loc[selected, "selected_component"] = winner
        else: winner, correlation = None, None
        phasing_qc.append({"region":v.name, "component":winner, "correlation":correlation, "status":"oriented" if correlation is not None and abs(correlation) >= c["MIN_PHASING_CORRELATION"] else "ambiguous"})
    prefix = Path(c["run_dir"])/f"03_features/compartments/{scope}"
    qcp = prefix.with_suffix(".phasing.json"); atomic_json(qcp, phasing_qc)
    return [save(vectors.drop(columns="value"), prefix.with_suffix(".tsv")), save(eigenvalues, prefix.with_suffix(".eigenvalues.tsv")), qcp, track(vectors,"oriented_score",prefix.with_suffix(".bedgraph"))]


def insulation(c, scope):
    import cooltools
    clr, bins, mask, labels, views, values = context(c, scope, c["INSULATION_RESOLUTION"])
    windows = [int(x) for x in c["INSULATION_WINDOWS_BP"].split(",")]
    table = cooltools.insulation(clr, windows, view_df=tool_views(clr,views), clr_weight_name="weight", ignore_diags=c["INSULATION_IGNORE_DIAGS"], min_frac_valid_pixels=c["INSULATION_MIN_VALID_PIXEL_FRACTION"], threshold=c["INSULATION_BOUNDARY_THRESHOLD"], append_raw_scores=True, chunksize=c["PIXEL_CHUNKSIZE"], nproc=c["THREADS_INSULATION"])
    table = registered_table(c, c["INSULATION_RESOLUTION"], table)
    for window in windows:
        width=window//c["INSULATION_RESOLUTION"]
        theoretical=sum(a+b>=c["INSULATION_IGNORE_DIAGS"] for a in range(width) for b in range(width))
        fraction=table[f"n_valid_pixels_{window}"]/theoretical if theoretical else np.zeros(len(table))
        table[f"valid_pixel_fraction_{window}"]=fraction
        insufficient=fraction<c["INSULATION_MIN_VALID_PIXEL_FRACTION"]
        table.loc[insufficient,[f"log2_insulation_score_{window}",f"boundary_strength_{window}"]]=np.nan
        table.loc[insufficient,f"is_boundary_{window}"]=False
    if "region" not in table:
        table_labels = view_labels(table, views)
        table["region"] = [views.iloc[label]["name"] if label >= 0 else pd.NA for label in table_labels]
    table["inference_status"] = "per_sample_measurement"
    outputs = [save(table, Path(c["run_dir"])/f"03_features/insulation/{scope}.tsv")]
    for window in windows:
        outputs.append(track(table,f"log2_insulation_score_{window}",Path(c["run_dir"])/f"03_features/insulation/{scope}.{window}.bedgraph"))
        boundary = table.loc[table[f"is_boundary_{window}"].fillna(False), ["chrom", "start", "end", "region", f"boundary_strength_{window}"]]
        outputs += [save(boundary, Path(c["run_dir"])/f"03_features/boundaries/{scope}.{window}.tsv")]
    return outputs


def loops(c, scope):
    outputs = []
    for res in loop_resolutions(c):
        clr, bins, mask, labels, views, values = context(c, scope, res)
        if c["LOOP_CALLER"] == "mustache":
            from .fine_loops import mustache_calls
            calls = mustache_calls(c,clr,bins,mask,labels,views,values,res)
        else:
            import cooltools
            expected_path = Path(c["run_dir"])/f"01_qc/ps/{scope}.{res}.tsv"
            expected = pd.read_csv(expected_path, sep="\t") if expected_path.exists() else expected_table(clr, bins, mask, labels, views, values, c)
            calls = cooltools.dots(clr, expected, expected_value_col="balanced.avg", clr_weight_name="weight", view_df=tool_views(clr,views), max_loci_separation=c["LOOP_MAX_DISTANCE_BP"], nproc=c["THREADS_LOOPS"])
        resource_plans = calls.attrs.get('resource_plans', [])
        from .loop_profiles import mustache_geometry
        geometry = mustache_geometry(c, res) if c['LOOP_CALLER'] == 'mustache' else None
        if len(calls):
            calls = calls.loc[(calls.chrom1 == calls.chrom2) & (calls.start2-calls.start1 >= c["LOOP_MIN_DISTANCE_BP"]) & (calls.start2-calls.start1 <= c["LOOP_MAX_DISTANCE_BP"])]
            calls = canonical_loops(calls, c["REFERENCE_ID"], res)
        else:
            calls = pd.DataFrame(columns=LOOP_COLUMNS)
        calls["caller"] = c["LOOP_CALLER"]; calls["resolution"] = res
        calls["inference_status"] = "discovery_enrichment_only"
        outputs += [save(calls, Path(c["run_dir"])/f"03_features/loops/{scope}.{res}.bedpe.tsv")]
        provenance=Path(c["run_dir"])/f"03_features/loops/{scope}.{res}.json"
        atomic_json(provenance,{"sample_id":scope,"resolution":res,"matrix_uri":clr.uri,"caller":c["LOOP_CALLER"],
                    "software_inventory":inventory(),"condition_pool":scope.startswith("pool__"),
                    "normalization":"inherited_stored_weights; common-mask excluded pixels unavailable",
                    "min_distance_bp":c["LOOP_MIN_DISTANCE_BP"],"max_distance_bp":c["LOOP_MAX_DISTANCE_BP"],
                    "caller_intrinsic_excluded_diagonals":4 if c["LOOP_CALLER"]=="mustache" else None,
                    "effective_minimum_candidate_bin_distance_bp":geometry['effective_best_case_min_distance_bp'] if geometry else c["LOOP_MIN_DISTANCE_BP"],
                    "geometry": geometry, "resource_plans": resource_plans,
                    "parameters":{k:v for k,v in c.items() if k.startswith("MUSTACHE_")},"profile":c["ANALYSIS_PROFILE"],
                    "methods_status":"adapted; not an exact Hsieh reproduction"})
        outputs.append(provenance)
    return outputs


def consensus(c, scope):
    root = Path(c["run_dir"])
    outputs = []
    sources=[]
    if c["RUN_LOOPS"]:
        for res in loop_resolutions(c):
            frames={sid:pd.read_csv(root/f"03_features/loops/{sid}.{res}.bedpe.tsv",sep="\t").assign(caller=c["LOOP_CALLER"]) for sid in stage_scopes(c,"loops")}
            sources.extend((sid,res,frame) for sid,frame in frames.items())
            outputs.append(save(consensus_loops(frames,c["REFERENCE_ID"],res,c["LOOP_CONSENSUS_TOLERANCE_BINS"]),root/f"03_features/loops/consensus.{res}.tsv"))
    if c.get("RUN_CHROMOSIGHT_DISCOVERY"):
        from .resolutions import chromosight_discovery_resolutions
        for res in chromosight_discovery_resolutions(c):
            frames = {sid: pd.read_csv(root/f"03_features/loops/{sid}.{res}.chromosight.tsv", sep="\t") for sid in stage_scopes(c, "chromosight_discovery")}
            sources.extend((sid, res, frame) for sid, frame in frames.items())
    from .fine_loops import cross_resolution_catalogue, map_candidates
    catalogue,members=cross_resolution_catalogue(c,sources)
    outputs += [save(catalogue,root/"03_features/loops/catalogue.tsv"),save(members,root/"03_features/loops/source_membership.tsv")]
    loops_frame,_=map_candidates(catalogue,c["REFERENCE_ID"],c["LOOP_RESOLUTION"])
    outputs += [save(loops_frame, root/"03_features/loops/consensus.tsv")]
    domains = []
    if c["RUN_INSULATION"]:
        with np.load(root/f"00_metadata/masks/{c['INSULATION_RESOLUTION']}.npz") as archive:
            common_mask = archive["mask"].copy()
        bins = pd.read_csv(root/f"00_metadata/masks/{c['INSULATION_RESOLUTION']}.bins.tsv", sep="\t")
        domain_index = masked_interval_index(bins, common_mask)
        for window in map(int, c["INSULATION_WINDOWS_BP"].split(",")):
            boundaries = {sid:pd.read_csv(root/f"03_features/boundaries/{sid}.{window}.tsv",sep="\t") for sid in stage_scopes(c,"insulation")}
            merged = boundary_consensus(boundaries, c["REFERENCE_ID"], c["INSULATION_RESOLUTION"], window)
            outputs += [save(merged, root/f"03_features/boundaries/consensus.{window}.tsv")]
            for (chrom, region), group in merged.groupby(["chrom", "region"]):
                positions = sorted(group.start)
                for start, end in zip(positions, positions[1:]):
                    if not interval_usable(domain_index, chrom, start, end): continue
                    domains.append({"chrom":chrom, "start":start, "end":end, "region":region, "window_bp":window,
                                    "feature_id":feature_id("domain",c["REFERENCE_ID"],(chrom,start,end,window),c["INSULATION_RESOLUTION"]), "inference_status":"descriptive"})
    outputs += [save(pd.DataFrame(domains, columns=["chrom","start","end","region","window_bp","feature_id","inference_status"]), root/"03_features/boundary_derived_domains/consensus.tsv")]
    return outputs


def masked_interval_index(bins, mask):
    result = {}
    for chrom, frame in bins.groupby("chrom", sort=False, observed=True):
        frame = frame.sort_values("start")
        starts, ends = frame.start.to_numpy(), frame.end.to_numpy()
        require((np.diff(ends) >= 0).all(), "Registered bin ends must be ordered")
        result[str(chrom)] = (starts, ends, np.r_[0, np.cumsum(~mask[frame.index])])
    return result


def interval_usable(index, chrom, start, end):
    if chrom not in index:
        return False
    starts, ends, bad = index[chrom]
    a = int(np.searchsorted(starts, start, side="left"))
    b = int(np.searchsorted(ends, end, side="right"))
    return b > a and bad[b] == bad[a]


def quantify_features(c, scope):
    root = Path(c["run_dir"])
    samples = validate_metadata(c)["samples"]
    discovered_path = root/"03_features/loops/catalogue.tsv"
    discovered = pd.read_csv(discovered_path,sep="\t") if has_loop_discovery(c) and discovered_path.exists() else pd.DataFrame(columns=LOOP_COLUMNS)
    independent, registry = external_features(c,samples)
    outputs=[]
    from .fine_loops import map_candidates
    for resolution in quantification_resolutions(c):
        measurements = {}
        for name, universe in (("independent_loops", independent), ("discovered_loops", discovered)):
            measurement, aliases = map_candidates(universe.reset_index(drop=True), c["REFERENCE_ID"], resolution)
            measurements[name] = (measurement, aliases)
        matrices = {name: frame[LOOP_COLUMNS].copy() for name, (frame, _) in measurements.items()}
        exposure_rows = []
        exclusions_by_name = {name: [] for name in measurements}
        from .quantification import quantify_observations
        sample_results = quantify_observations(c, samples, resolution, measurements)
        for sample, results in zip(samples, sample_results):
            for name, (counts, exposure, reasons) in results.items():
                matrices[name][sample["sample_id"]] = counts
                measurement = measurements[name][0]
                exclusions_by_name[name].extend({"sample_id": sample["sample_id"], "feature_id": measurement.loc[i, "feature_id"], "reason": reason} for i, reason in reasons.items())
            exposure_rows.append({"sample_id": sample["sample_id"], "eligible_cis_exposure": exposure, "resolution": resolution})
        for name, (measurement, aliases) in measurements.items():
            matrix = matrices[name]
            exclusions = exclusions_by_name[name]
            directory=root/"04_quantification"/name/str(resolution)
            tables={"counts.tsv":matrix,"exposures.tsv":pd.DataFrame(exposure_rows),"exclusions.tsv":pd.DataFrame(exclusions,columns=["sample_id","feature_id","reason"]),"candidate_mapping.tsv":aliases}
            for filename,table in tables.items():
                outputs.append(save(table,directory/filename))
                if resolution==c["LOOP_RESOLUTION"]: outputs.append(save(table,directory.parent/filename))
            info={"inference_status":"independent_candidates" if name=="independent_loops" and registry else "descriptive_same_data", "registry":registry if name=="independent_loops" else None, "feature_hash":digest(measurement[LOOP_COLUMNS].to_dict("records")),"resolution":resolution,"collapsed_candidates_are_one_measurement":True}
            provenance=directory/"features.json";atomic_json(provenance,info);outputs.append(provenance)
            if resolution==c["LOOP_RESOLUTION"]:
                provenance=directory.parent/"features.json";atomic_json(provenance,info);outputs.append(provenance)
    return outputs


def expected_lookup(expected):
    """Build distance arrays once for an immutable sample/resolution context.

    This is an in-memory lookup, never a cache reused across runs or masks.
    """
    if isinstance(expected, dict):
        return expected
    result = {}
    for (left, right), frame in expected.groupby(["region1", "region2"], sort=False):
        if left != right:
            continue
        distances = frame.dist.to_numpy(dtype=int)
        require((distances >= 0).all(), "Expected contact distances must be nonnegative")
        values = np.full(int(distances.max())+1 if len(distances) else 0, np.nan)
        values[distances] = frame["balanced.avg"].to_numpy()
        result[left] = values
    return result


def oe_window(clr, region1, region2, values, labels, views, expected, c):
    lo,hi=clr.extent(region1); lo2,hi2=clr.extent(region2)
    raw=bounded_raw(clr,region1,region2,c).astype(float)
    result=raw*values[lo:hi,None]*values[None,lo2:hi2]
    i=np.arange(lo,hi)[:,None]; j=np.arange(lo2,hi2)[None,:]
    distance=np.abs(j-i)
    denominator=np.full(result.shape,np.nan)
    lookup_by_view = expected_lookup(expected)
    for vid in np.unique(labels[lo:hi]):
        if vid<0:continue
        name=views.iloc[vid]["name"]
        lookup = lookup_by_view.get(name, np.empty(0))
        same=(labels[lo:hi,None]==vid)&(labels[None,lo2:hi2]==vid)
        valid = same & (distance < len(lookup))
        denominator[valid]=lookup[distance[valid]]
    good=np.isfinite(denominator)&(denominator>0)&np.isfinite(result)
    return np.divide(result,denominator,out=np.full(result.shape,np.nan),where=good)


def saddle_partition(c, samples):
    root=Path(c["run_dir"])
    resolution=c["COMPARTMENT_RESOLUTION"]
    base=pd.read_csv(root/f"00_metadata/masks/{resolution}.bins.tsv",sep="\t")
    archive=np.load(root/f"00_metadata/masks/{resolution}.npz")
    labels=archive["labels"]
    tracks=[]
    if c["SADDLE_PARTITION_SOURCE"]=="external":
        track=pd.read_csv(c["SADDLE_REFERENCE_TRACK"],sep="\t",header=None,names=["chrom","start","end","oriented_score"])
        joined=base.merge(track,on=["chrom","start","end"],how="left",validate="one_to_one")
        tracks=[joined.oriented_score.to_numpy()]
    else:
        if c["SADDLE_PARTITION_SOURCE"]=="all_sample_pool":
            ids=["pool__"+p["pool_id"] for p in read_json(root/"02_pools/registry.json")["pools"] if p["pool_scope"]=="all_samples"]
        else:
            ids=[s["sample_id"] for s in samples]
            require(len(ids)==len(set(ids)), "Saddle observation IDs must be unique")
            manifest_order={s["sample_id"]:i for i,s in enumerate(validate_metadata(c)["samples"])}
            require(set(ids)<=set(manifest_order), "Saddle observation IDs must belong to the current manifest")
            # Floating-point rank means can differ at ties when tracks are added
            # in another order. Preserve normal manifest-order numerics for every
            # caller ordering before the existing coordinate tie-break.
            ids.sort(key=manifest_order.__getitem__)
        require(ids,"Requested all-sample pool is unavailable")
        for sample in ids:
            table=pd.read_csv(root/f"03_features/compartments/{sample}.tsv",sep="\t")
            joined=base.merge(table[["chrom","start","end","oriented_score"]],on=["chrom","start","end"],how="left",validate="one_to_one")
            tracks.append(joined.oriented_score.to_numpy())
    valid=archive["mask"]&np.isfinite(np.stack(tracks)).all(axis=0)
    assignments=np.full(len(base),-1,dtype=int)
    score=np.full(len(base),np.nan)
    for vid in np.unique(labels[valid]):
        ids=np.flatnonzero(valid&(labels==vid))
        if len(ids)<50:continue
        ranks=np.stack([pd.Series(t[ids]).rank(method="average",pct=True).to_numpy() for t in tracks])
        score[ids]=ranks.mean(axis=0)
        ordered=ids[np.lexsort((ids,score[ids]))]
        assignments[ordered]=np.minimum(np.arange(len(ids))*50//len(ids),49)
    base["partition"]=assignments; base["reference_rank"]=score
    return base,labels


def _aggregates_single(c,scope):
    modules=aggregate_modules(c)
    root=Path(c["run_dir"])
    apa_directory=root/'05_aggregates/apa'
    if c.get('_apa_resolution'):apa_directory=apa_directory/str(c['_apa_resolution'])
    samples=validate_metadata(c)["samples"]
    loops_path=root/"03_features/loops/consensus.tsv"
    discovered=pd.read_csv(loops_path,sep="\t") if "apa" in modules and has_loop_discovery(c) and loops_path.exists() else pd.DataFrame(columns=LOOP_COLUMNS)
    independent,_=external_features(c,samples) if "apa" in modules and not c.get("_apa_ignore_external_registry") else (pd.DataFrame(columns=LOOP_COLUMNS), None)
    universe=c.get("_apa_candidates",independent if len(independent) else discovered)
    domains_path=root/"03_features/boundary_derived_domains/consensus.tsv"
    domains=pd.read_csv(domains_path,sep="\t") if "ada" in modules and c["RUN_INSULATION"] and domains_path.exists() else pd.DataFrame()
    if len(domains):domains=domains.loc[domains.window_bp==c["INSULATION_PRIMARY_WINDOW_BP"]]
    outputs=[]; summaries=[]; domain_rows=[]
    selected_samples = [s for s in samples if not c.get("_aggregate_sample_id") or s["sample_id"] == c["_aggregate_sample_id"]]
    unit_dir = Path(c["_aggregate_unit_dir"]) if c.get("_aggregate_unit_dir") else None
    from .aggregation import atomic_npz
    from .expected_cache import expected_cached
    for sample in selected_samples if modules & {"apa", "ada"} else []:
        sid=sample["sample_id"];res=c["LOOP_RESOLUTION"]
        clr,bins,mask,labels,views,values=context(c,sid,res)
        expected=expected_lookup(expected_cached(c,sid,res,clr,bins,mask,labels,views,values))
        if "apa" in modules:
            flank=c["APA_FLANK_BP"]//res; shape=2*flank+1
            total=np.zeros((shape,shape));coverage=np.zeros((shape,shape),dtype=int)
            for feature in universe.itertuples():
                region1=(feature.chrom1,feature.start1-c["APA_FLANK_BP"],feature.end1+c["APA_FLANK_BP"])
                region2=(feature.chrom2,feature.start2-c["APA_FLANK_BP"],feature.end2+c["APA_FLANK_BP"])
                if region1[1]<0 or region2[1]<0 or region1[2]>clr.chromsizes[feature.chrom1] or region2[2]>clr.chromsizes[feature.chrom2]:continue
                window=oe_window(clr,region1,region2,values,labels,views,expected,c)
                if window.shape!=total.shape:continue
                good=np.isfinite(window)
                total[good]+=window[good];coverage[good]+=1
            average=np.divide(total,coverage,out=np.full(total.shape,np.nan),where=coverage>0)
            path=apa_directory/f"{sid}.npz";path.parent.mkdir(parents=True,exist_ok=True)
            atomic_npz(path,mean=average,sum=total,valid_feature_count=coverage)
            outputs.append(path)
            summaries.append({"sample_id":sid,"resolution":res,"requested_features":len(universe),"central_oe":average[flank,flank],"central_valid_features":coverage[flank,flank],"central_3x3_oe":float(np.nanmean(average[flank-1:flank+2,flank-1:flank+2])) if np.isfinite(average).any() else np.nan,"inference_status":"descriptive"})
        if "ada" in modules:
            ada_total=np.zeros((50,50));ada_count=np.zeros((50,50),dtype=int)
            for feature in domains.itertuples():
                window=oe_window(clr,(feature.chrom,feature.start,feature.end),(feature.chrom,feature.start,feature.end),values,labels,views,expected,c)
                upper=np.triu(np.ones(window.shape,dtype=bool),k=c["INSULATION_IGNORE_DIAGS"])
                finite=upper&np.isfinite(window)
                denominator=int(finite.sum());total_oe=float(window[finite].sum())
                domain_rows.append({"feature_id":feature.feature_id,"sample_id":sid,"mean_intra_domain_oe":total_oe/denominator if denominator else np.nan,"oe_sum":total_oe,"eligible_pixel_count":denominator,"inference_status":"descriptive"})
                # Bin-average rescaling; zero and missing pixels remain distinguishable.
                cuts=np.linspace(0,len(window),51,dtype=int)
                for i,j in itertools.product(range(50),repeat=2):
                    block=window[cuts[i]:cuts[i+1],cuts[j]:cuts[j+1]]
                    if np.isfinite(block).any():ada_total[i,j]+=float(np.nanmean(block));ada_count[i,j]+=1
            path=root/f"05_aggregates/ada/{sid}.npz";path.parent.mkdir(parents=True,exist_ok=True)
            atomic_npz(path,mean=np.divide(ada_total,ada_count,out=np.full((50,50),np.nan),where=ada_count>0),sum=ada_total,valid_domain_count=ada_count)
            outputs.append(path)
    if "apa" in modules:
        outputs.append(save(pd.DataFrame(summaries),(unit_dir/"summary.tsv") if unit_dir else apa_directory/"summary.tsv"))
    if "ada" in modules:
        outputs.append(save(pd.DataFrame(domain_rows,columns=["feature_id","sample_id","mean_intra_domain_oe","oe_sum","eligible_pixel_count","inference_status"]),(unit_dir/"strength.tsv") if unit_dir else root/"04_quantification/domain_scores/strength.tsv"))
    if "saddle" in modules:
        partition,labels=c["_saddle_partition"] if c.get("_saddle_partition") is not None else saddle_partition(c,samples)
        require((partition.partition >= 0).any(),
                "Saddle has no eligible partition: each usable view needs at least 50 common, consistently oriented bins; inspect compartment phasing QC and the registered track/grid")
        if unit_dir is None: outputs.append(save(partition,root/"05_aggregates/saddle/partition.tsv"))
        info=root/"05_aggregates/saddle/partition.json"
        if unit_dir is None: atomic_json(info,{"source":c["SADDLE_PARTITION_SOURCE"],"groups":50,"partition_hash":digest(partition[["chrom","start","end","partition"]].to_dict("records")),"views":[{"view_index":int(vid),"usable_bins":int(((labels==vid)&(partition.partition.to_numpy()>=0)).sum()),"status":"partitioned" if ((labels==vid)&(partition.partition.to_numpy()>=0)).sum()>=50 else "fewer_than_50_usable_bins"} for vid in np.unique(labels[labels>=0])],"inference_status":"descriptive"})
        if unit_dir is None: outputs.append(info)
        saddle_rows=[]
        for sample in selected_samples:
            sid=sample["sample_id"];clr,bins,mask,labels,views,values=context(c,sid,c["COMPARTMENT_RESOLUTION"])
            ex=expected_lookup(expected_cached(c,sid,c["COMPARTMENT_RESOLUTION"],clr,bins,mask,labels,views,values))
            assignments=np.full(len(bins),-1,dtype=int)
            assignments[native_indices(clr,partition[["chrom","start","end"]])]=partition.partition.to_numpy()
            sums=np.zeros((50,50));counts=np.zeros((50,50),dtype=np.int64)
            # Dense per-view only after a hard bound; no genome-wide dense matrix.
            for vid,view in enumerate(views.itertuples()):
                ids=np.flatnonzero(labels==vid)
                region=(view.chrom,view.start,view.end)
                matrix=oe_window(clr,region,region,values,labels,views,ex,c)
                groups=assignments[ids]
                for a,b in itertools.product(range(50),repeat=2):
                    if not (groups==a).any() or not (groups==b).any():continue
                    block=matrix[np.ix_(groups==a,groups==b)]
                    good=np.isfinite(block)
                    sums[a,b]+=np.nansum(block);counts[a,b]+=good.sum()
            average=np.divide(sums,counts,out=np.full((50,50),np.nan),where=counts>0)
            path=root/f"05_aggregates/saddle/{sid}.npz";path.parent.mkdir(parents=True,exist_ok=True)
            atomic_npz(path,mean=average,sum=sums,valid_pixel_count=counts);outputs.append(path)
            def corner(x,y):
                n=counts[x,y].sum();return sums[x,y].sum()/n if n else np.nan
            aa=corner(slice(40,50),slice(40,50));bb=corner(slice(0,10),slice(0,10));ab=corner(slice(0,10),slice(40,50))
            saddle_rows.append({"sample_id":sid,"AA":aa,"BB":bb,"AB":ab,"segregation_ratio":(aa+bb)/(2*ab) if np.isfinite(ab) and ab>0 else np.nan,"inference_status":"descriptive"})
        outputs.append(save(pd.DataFrame(saddle_rows),(unit_dir/"summary.tsv") if unit_dir else root/"05_aggregates/saddle/summary.tsv"))
    return outputs


def aggregates(c,scope):
    """Independent ADA/saddle branches and a complete APA path per loop grid."""
    modules=aggregate_modules(c);root=Path(c["run_dir"]);outputs=[]
    for module in sorted(modules):
        resolutions=quantification_resolutions(c) if module=="apa" else [c["ADA_RESOLUTION"]] if module=="ada" else [c["COMPARTMENT_RESOLUTION"]]
        for res in resolutions:
            config=dict(c,AGGREGATE_MODULES=module,LOOP_RESOLUTION=res)
            if module=="apa":
                config['_apa_resolution']=res
                # Fixed cross-resolution universe, mapped once per measured pixel.
                catalogue=root/"03_features/loops/catalogue.tsv"
                from .fine_loops import map_candidates
                if has_loop_discovery(c) and catalogue.exists():
                    mapped,_=map_candidates(pd.read_csv(catalogue,sep="\t"),c["REFERENCE_ID"],res)
                    config["_apa_candidates"]=mapped
                independent,_=external_features(c,validate_metadata(c)["samples"])
                if len(independent): config["_apa_candidates"]=map_candidates(independent,c["REFERENCE_ID"],res)[0]
                config["_apa_ignore_external_registry"]=True
            from .aggregation import run_units, atomic_npz
            products=run_units(config,module,res)
            if module=="apa":
                for path in products:
                    outputs.append(path)
                    if res==c["LOOP_RESOLUTION"]:
                        alias=path.parent.parent/path.name;shutil.copyfile(path,alias);outputs.append(alias)
            else: outputs.extend(products)
            if module == "ada":
                domain_path = root/"03_features/boundary_derived_domains/consensus.tsv"
                requested = pd.read_csv(domain_path, sep="\t")
                requested = requested.loc[requested.window_bp == c["INSULATION_PRIMARY_WINDOW_BP"]].copy()
                outputs.append(save(requested, root/"05_aggregates/ada/requested_domains.tsv"))
                contribution_path = root/"05_aggregates/ada/contributions.tsv"
                scores = pd.read_csv(root/"04_quantification/domain_scores/strength.tsv", sep="\t")
                scores['status'] = np.where(scores.eligible_pixel_count > 0, 'eligible_scalar_pixels', 'no_eligible_scalar_pixels')
                outputs.append(save(scores, contribution_path))
            location=root/"05_aggregates"/module/(str(res) if module=="apa" else "")
            samples=validate_metadata(c)["samples"]
            for condition in sorted({s["condition"] for s in samples}):
                members=[s for s in samples if s["condition"]==condition]
                arrays=[np.load(location/f"{s['sample_id']}.npz")["mean"] for s in members]
                stack=np.stack(arrays);valid=np.isfinite(stack);count=valid.sum(axis=0)
                average=np.divide(np.nansum(stack,axis=0),count,out=np.full(stack.shape[1:],np.nan),where=count>0)
                condition_id=condition if IDENTIFIER.fullmatch(condition) else 'condition_'+digest(condition)[:12]
                path=location/("condition__"+condition_id+".npz")
                atomic_npz(path,mean=average,valid_replicate_count=count);outputs.append(path)
            metadata_path=location/"aggregation.json"
            atomic_json(metadata_path,{"module":module,"resolution":res,"condition_summary":"equal biological-replicate weighting; missing pixels excluded, never zero-filled","condition_members":{condition:[s["sample_id"] for s in samples if s["condition"]==condition] for condition in sorted({s["condition"] for s in samples})},"inference_status":"descriptive"});outputs.append(metadata_path)
    return outputs


def hicrep(c, scope, derived_root=None):
    from scipy.sparse import diags
    from hicrep.hicrep import sccByDiag, varVstran
    from hicrep.utils import meanFilterSparse, trimDiags, upperDiagCsr
    samples=validate_metadata(c)["samples"]
    rows=[]
    for left,right in itertools.combinations(samples,2):
        clr,bins,mask,labels,views,_=context(c,left["sample_id"],c["HICREP_RESOLUTION"])
        other,*_=context(c,right["sample_id"],c["HICREP_RESOLUTION"])
        if derived_root is not None:
            clr=cooler.Cooler(str(Path(derived_root)/left["sample_id"]/f"{c['HICREP_RESOLUTION']}.cool"))
            other=cooler.Cooler(str(Path(derived_root)/right["sample_id"]/f"{c['HICREP_RESOLUTION']}.cool"))
            base=Path(c["run_dir"])/f"00_metadata/masks/{c['HICREP_RESOLUTION']}"
            common_bins=pd.read_csv(str(base)+".bins.tsv",sep="\t")
            with np.load(str(base)+".npz") as archive:
                _,mask,_,_=native_domain(clr,common_bins,archive["mask"],archive["labels"],np.ones(len(common_bins)))
        for vid,view in enumerate(views.itertuples()):
            lo,hi=clr.extent((view.chrom,view.start,view.end))
            local=mask[lo:hi]
            matrices=[]
            for source in (clr,other):
                matrix=source.matrix(balance=False,sparse=True).fetch((view.chrom,view.start,view.end)).tocoo()
                matrix=(diags(local.astype(float))@matrix@diags(local.astype(float))).tocoo()
                n_diags=min(matrix.shape[0],c["HICREP_MAX_DISTANCE_BP"]//clr.binsize+1)
                matrix=trimDiags(matrix,n_diags,False).astype(float)
                total=float(matrix.sum())
                if total>0:matrix=matrix/total
                matrix=meanFilterSparse(matrix,c["HICREP_SMOOTHING_BINS"]).tocoo()
                # Exclude windows whose smoothing footprint touches an excluded bin.
                valid=local.copy()
                for offset in range(1,c["HICREP_SMOOTHING_BINS"]+1):
                    valid[offset:] &= local[:-offset];valid[:-offset] &= local[offset:]
                matrix=(diags(valid.astype(float))@matrix@diags(valid.astype(float))).tocoo()
                matrices.append(trimDiags(matrix,n_diags,False))
            m1,m2=matrices
            score=float(sccByDiag(m1,m2,n_diags)) if m1.nnz and m2.nnz else np.nan
            d1=upperDiagCsr(m1,n_diags);d2=upperDiagCsr(m2,n_diags)
            n=(d1+d2).getnnz(axis=1)
            with np.errstate(divide="ignore",invalid="ignore"):
                weight=float(np.nan_to_num(n*varVstran(n),nan=0,posinf=0,neginf=0).sum())
            rows.append({"sample1":left["sample_id"],"sample2":right["sample_id"],"region":view.name,"scc":score,"stratum_weight":weight,"excluded_bin_count":int((~local).sum()),"smoothing_half_width_bins":c["HICREP_SMOOTHING_BINS"],"status":"valid" if np.isfinite(score) else "not_applicable_sparse"})
    frame=pd.DataFrame(rows,columns=["sample1","sample2","region","scc","stratum_weight","excluded_bin_count","smoothing_half_width_bins","status"])
    summary=[]
    for (left,right),group in frame.groupby(["sample1","sample2"]):
        valid=np.isfinite(group.scc)&(group.stratum_weight>0)
        score=np.average(group.loc[valid,"scc"],weights=group.loc[valid,"stratum_weight"]) if valid.any() else np.nan
        summary.append({"sample1":left,"sample2":right,"weighted_scc":score,"valid_views":int(valid.sum())})
    root=Path(c["run_dir"]) if derived_root is None else Path(derived_root)
    return [save(frame,root/"01_qc/hicrep/by_view.tsv"),save(pd.DataFrame(summary,columns=["sample1","sample2","weighted_scc","valid_views"]),root/"01_qc/hicrep/summary.tsv")]


def balance_new_map(c,path,blacklist=None,cis_only=False):
    require(c["BALANCE_NEW_DERIVED_MAPS"],"New derived-map balancing must be explicitly enabled")
    path=Path(str(path).split("::")[0]).resolve()
    root=Path(c["run_dir"]).resolve()
    require(path.is_relative_to(root/"02_pools") or path.is_relative_to(root/"01_qc/depth_sensitivity"),"Balancing is restricted to task-owned new derived maps")
    inputs=read_json(root/"00_metadata/inputs.json")
    require(str(path) not in {v["path"] for v in inputs.values()},"Imported maps must never be balanced")
    clr=cooler.Cooler(str(path))
    for sample in validate_metadata(c)["samples"]:
        source=cooler.Cooler(matrix_uri(sample_matrix(c,sample,clr.binsize),clr.binsize))
        if same_raw_counts(clr,source,c["PIXEL_CHUNKSIZE"]):
            raise ValueError("Derived map has unchanged imported counts; reuse its stored weights instead of balancing a copy")
    from multiprocessing import Pool
    with Pool(c["THREADS_BALANCE"]) as workers:
        weight,stats=cooler.balance_cooler(clr,store=True,blacklist=blacklist,cis_only=cis_only,chunksize=c["PIXEL_CHUNKSIZE"],ignore_diags=c["INSULATION_IGNORE_DIAGS"],map=workers.map)
    require(bool(np.asarray(stats["converged"]).all()),"New derived-map balancing did not converge")
    return stats


def same_raw_counts(left,right,chunk_size):
    if left.binsize!=right.binsize:
        return False
    if left.info.get("storage-mode")==right.info.get("storage-mode")=="symmetric-upper" and left.info["sum"]!=right.info["sum"]:return False
    axes=[]
    for clr in (left,right):
        table=clr.bins()[:][["chrom","start","end"]].copy();table["chrom"]=table.chrom.astype(str);axes.append(table.to_numpy())
    if not np.array_equal(*axes):return False
    def entries(clr):
        for chunk in pixel_chunks(clr,chunk_size):
            for row in chunk[["bin1_id","bin2_id","count"]].itertuples(index=False,name=None):
                if row[2]:yield row
    sentinel=object()
    return all(a==b for a,b in itertools.zip_longest(entries(left),entries(right),fillvalue=sentinel))


def validate_pool_sum(raw, members, chunk_size):
    """Compare every pooled pixel to its member sum without constructing another map."""
    axes=raw.bins()[:][["chrom","start","end"]].astype({"chrom":str})
    require(all(clr.bins()[:][["chrom","start","end"]].astype({"chrom":str}).equals(axes) for clr in members),
            "Imported pool raw-sum verification requires identical native member/pool grids; use explicit POOL_MODE=create for differing scopes")
    def entries(clr):
        for chunk in pixel_chunks(clr,chunk_size):
            for i,j,y in chunk[["bin1_id","bin2_id","count"]].itertuples(index=False,name=None):
                if y:yield (int(i),int(j)),int(y)
    combined=heapq.merge(*(entries(clr) for clr in members),key=lambda x:x[0])
    wanted=((key,sum(value for _,value in group)) for key,group in itertools.groupby(combined,key=lambda x:x[0]))
    sentinel=object()
    for actual,expected in itertools.zip_longest(entries(raw),wanted,fillvalue=sentinel):
        require(actual==expected,"Imported pool pixel counts differ from exact member raw sums; request corrected upstream pooling")


def pools(c,scope):
    root=Path(c["run_dir"]);samples=validate_metadata(c)["samples"]
    registered=[];outputs=[]
    if c["POOLS_MANIFEST"]:
        for row in read_table(c["POOLS_MANIFEST"],("pool_id","pool_scope","condition","member_sample_ids","reference_id","matrix_path","provenance_path")):
            require(IDENTIFIER.fullmatch(row["pool_id"]),"Unsafe pool_id")
            row["matrix_path"]=resolve(row["matrix_path"],Path(c["POOLS_MANIFEST"]).parent)
            row["provenance_path"]=resolve(row["provenance_path"],Path(c["POOLS_MANIFEST"]).parent)
            members=row["member_sample_ids"].split(";")
            require(len(set(members))==len(members) and set(members)<={s["sample_id"] for s in samples},"Pool membership invalid")
            require(row["reference_id"]==c["REFERENCE_ID"] and row["pool_scope"] in ("condition","all_samples"),"Pool reference/scope mismatch")
            expected_members=[s["sample_id"] for s in samples if row["pool_scope"]=="all_samples" or s["condition"]==row["condition"]]
            require(set(members)==set(expected_members),"Pool membership must exactly match its declared scope")
            prov=read_json(row["provenance_path"])
            require(prov.get("reference_id")==c["REFERENCE_ID"] and prov.get("count_units")=="unique_contacts","Pool reference/count provenance mismatch")
            require(set(prov["member_sample_ids"])==set(members) and prov.get("summation")=="raw_integer_sum", "Imported pool requires raw summation provenance")
            input_records=read_json(root/"00_metadata/inputs.json")
            require(all(prov["member_hashes"].get(s["sample_id"])==input_records[s["matrix_path"].split("::")[0]]["sha256"] for s in samples if s["sample_id"] in members),"Pool member fingerprint mismatch")
            source_prov=read_json(next(s["provenance_path"] for s in samples if s["sample_id"]==expected_members[0]))
            require(prov.get("filters")==source_prov["filters"],"Pool upstream filters differ from members")
            row["uris"]={}
            for res in resolution_set(c):
                raw=cooler.Cooler(matrix_uri(row["matrix_path"],res))
                bins=pd.read_csv(root/f"00_metadata/masks/{res}.bins.tsv",sep="\t")
                native_indices(raw,bins)
                weights(raw,c,prov,res)
                total=sum(int(p["count"].sum()) for p in pixel_chunks(raw,c["PIXEL_CHUNKSIZE"]))
                wanted=sum(unordered_total(cooler.Cooler(matrix_uri(sample_matrix(c,s,res),res)),c["PIXEL_CHUNKSIZE"]) for s in samples if s["sample_id"] in members)
                require(total==wanted,"Imported pool raw total differs from member sum")
                validate_pool_sum(raw,[cooler.Cooler(matrix_uri(sample_matrix(c,s,res),res)) for s in samples if s["sample_id"] in members],c["PIXEL_CHUNKSIZE"])
                row["uris"][str(res)]=raw.uri
            row["task_owned"]=False;registered.append(row)
    if c["POOL_MODE"]=="create":
        scopes=[("all_samples","",samples)]+[("condition",condition,[s for s in samples if s["condition"]==condition]) for condition in sorted({s["condition"] for s in samples})]
        for kind,condition,members in scopes:
            if any(p["pool_scope"]==kind and p["condition"]==condition for p in registered):continue
            pool_id="all_samples" if kind=="all_samples" else "condition_"+condition
            require(all(ch.isalnum() or ch in "_.-" for ch in pool_id),"Pool condition requires a safe identifier")
            if len(members)==1:
                sample=members[0]
                registered.append({"pool_id":pool_id,"pool_scope":kind,"condition":condition,"member_sample_ids":sample["sample_id"],"provenance_path":sample["provenance_path"],"task_owned":False,"uris":{str(res):matrix_uri(sample_matrix(c,sample,res),res) for res in resolution_set(c)}})
                continue
            require(c["BALANCE_NEW_DERIVED_MAPS"],"Creating pools for balanced discovery requires explicit new-map balancing")
            row={"pool_id":pool_id,"pool_scope":kind,"condition":condition,"member_sample_ids":";".join(s["sample_id"] for s in members),"task_owned":True,"uris":{}}
            input_records=read_json(root/"00_metadata/inputs.json")
            source_prov=read_json(members[0]["provenance_path"])
            prov={"reference_id":c["REFERENCE_ID"],"count_units":"unique_contacts","filters":source_prov["filters"],"balancing":{},"member_sample_ids":[s["sample_id"] for s in members],"member_hashes":{s["sample_id"]:input_records[s["matrix_path"].split("::")[0]]["sha256"] for s in members},"dedup_scope":"upstream_biological_members_already_deduplicated","library_ids":[lib for s in members for lib in read_json(s["provenance_path"])["library_ids"]],"source_ids":[source for s in members for source in read_json(s["provenance_path"])["source_ids"]],"summation":"raw_integer_sum"}
            for res in resolution_set(c):
                path=root/f"02_pools/{pool_id}/{res}.cool";path.parent.mkdir(parents=True,exist_ok=True)
                uris=[matrix_uri(sample_matrix(c,s,res),res) for s in members]
                pending=path.with_name("."+path.name+".pending.cool")
                bins=pd.read_csv(root/f"00_metadata/masks/{res}.bins.tsv",sep="\t")
                expected_sum=[0]
                def summed_chunks():
                    for uri in uris:
                        for chunk in aligned_pixel_chunks(cooler.Cooler(uri),bins,c["PIXEL_CHUNKSIZE"]):
                            expected_sum[0]+=int(chunk["count"].sum())
                            yield chunk
                scratch=Path(c["SCRATCH_ROOT"]);scratch.mkdir(parents=True,exist_ok=True)
                cooler.create_cooler(str(pending),bins,summed_chunks(),ordered=False,
                                     mergebuf=c["PIXEL_CHUNKSIZE"],temp_dir=str(scratch),
                                     dtypes={"count":"int64"},assembly=c["REFERENCE_ID"])
                require(cooler.Cooler(str(pending)).info["sum"]==expected_sum[0],"New pool violates selected-coordinate raw sum conservation")
                mask=np.load(root/f"00_metadata/masks/{res}.npz")["mask"]
                balance_new_map(c,pending,blacklist=np.flatnonzero(~mask))
                os.replace(pending,path)
                prov["balancing"][str(res)]={"weight_column":"weight","convention":"multiplicative","converged":True,"recipe":{"derived":"raw_pool"}}
                row["uris"][str(res)]=str(path);outputs.append(path)
            path=root/f"02_pools/{pool_id}/provenance.json";atomic_json(path,prov);row["provenance_path"]=str(path);outputs.append(path);registered.append(row)
    require(registered,"Enabled pool stage produced no registered pools")
    require(len({r["pool_id"] for r in registered})==len(registered),"Duplicate pool IDs")
    registry=root/"02_pools/registry.json";atomic_json(registry,{"pools":registered});outputs.append(registry)
    return outputs


def run_r(c,script,args):
    script=str(Path(__file__).parents[2]/"scripts"/script)
    if c["ENVIRONMENT_MODE"]=="development":command=[c["R_SCRIPT"],"--vanilla",script,*map(str,args)]
    else:command=[c["CONDA_EXE"],"run","--no-capture-output","-p",c["DIFF_R_ENV_PREFIX"],"Rscript","--vanilla",script,*map(str,args)]
    subprocess.run(command,check=True)


def _differential_single(c,scope):
    root=Path(c["run_dir"]);metadata=validate_metadata(c);outputs=[]
    for comparison in metadata["comparisons"]:
        if c.get("_comparison_modules"):
            comparison = dict(comparison, modules=";".join(m for m in comparison["modules"].split(";") if m in c["_comparison_modules"]))
            if not comparison["modules"]: continue
        design=design_matrix(metadata["samples"],comparison,c)
        ids=design["sample_ids"]
        case=[s["sample_id"] for s in metadata["samples"] if s["condition"]==comparison["case_condition"]]
        control=[s["sample_id"] for s in metadata["samples"] if s["condition"]==comparison["reference_condition"]]
        directory=root/"06_comparisons"/comparison["comparison_id"]
        if c.get("_loop_measurement_resolution"):directory=directory/str(c["_loop_measurement_resolution"])
        directory.mkdir(parents=True,exist_ok=True)
        status=check_calibration(c,comparison,read_json(root/"00_metadata/environments.json"))
        parameters={"design":design,"validation_status":status,"inference":comparison["inference"],"fdr":c["FDR_THRESHOLD"],"min_total":c["LOOP_MIN_TOTAL_COUNT"],"min_mean_at_min_exposure":c["LOOP_MIN_MEAN_COUNT_AT_MIN_EXPOSURE"],"min_features":c["LOOP_MIN_FEATURES_FORMAL"]}
        params=directory/"parameters.json";atomic_json(params,parameters);outputs.append(params)
        for module in comparison["modules"].split(";"):
            output=directory/(module+".tsv")
            if module=="loops":
                kind="independent_loops" if comparison["inference"]!="descriptive" else "discovered_loops"
                measurement=str(c.get("_loop_measurement_resolution",c["LOOP_RESOLUTION"]))
                counts=pd.read_csv(root/f"04_quantification/{kind}/{measurement}/counts.tsv",sep="\t")
                if comparison["inference"]=="descriptive" and not len(counts):
                    kind="independent_loops"
                    counts=pd.read_csv(root/f"04_quantification/{kind}/{measurement}/counts.tsv",sep="\t")
                parameters.update(loop_universe=kind,measurement_resolution=int(measurement))
                atomic_json(params,parameters)
                exposures=pd.read_csv(root/f"04_quantification/{kind}/{measurement}/exposures.tsv",sep="\t").set_index("sample_id").eligible_cis_exposure
                model_input=directory/"loop_model_counts.tsv";save(counts[LOOP_COLUMNS+ids],model_input);outputs.append(model_input)
                exposure_input=directory/"loop_model_exposures.tsv";save(exposures.loc[ids].rename_axis("sample_id").reset_index(),exposure_input);outputs.append(exposure_input)
                if comparison["inference"]!="descriptive":run_r(c,"differential_loops.R",[model_input,exposure_input,params,output])
                else:
                    signal=counts[ids].div(exposures[ids],axis=1)
                    result=counts[LOOP_COLUMNS].copy();result["case_mean_allocation"]=signal[case].mean(axis=1);result["reference_mean_allocation"]=signal[control].mean(axis=1)
                    with np.errstate(divide="ignore",invalid="ignore"):result["log2_allocation_ratio"]=np.log2(result.case_mean_allocation/result.reference_mean_allocation)
                    result["PValue"]=np.nan;result["FDR"]=np.nan;result["inference_status"]="descriptive";save(result,output)
            elif module=="insulation":
                matrices=[]
                for sid in ids:
                    frame=pd.read_csv(root/f"03_features/insulation/{sid}.tsv",sep="\t")
                    col=f"log2_insulation_score_{c['INSULATION_PRIMARY_WINDOW_BP']}"
                    matrices.append(frame[["chrom","start","end",col]].rename(columns={col:sid}))
                merged=matrices[0]
                for matrix in matrices[1:]:merged=merged.merge(matrix,on=["chrom","start","end"],how="outer",validate="one_to_one")
                annotations=[]
                for sid in ids:
                    table=pd.read_csv(root/f"03_features/insulation/{sid}.tsv",sep="\t")
                    annotations.append(table[["chrom","start","end",f"is_boundary_{c['INSULATION_PRIMARY_WINDOW_BP']}"]].rename(columns={f"is_boundary_{c['INSULATION_PRIMARY_WINDOW_BP']}":"called"}))
                calls=pd.concat(annotations).groupby(["chrom","start","end"],as_index=False).called.sum().rename(columns={"called":"boundary_support"})
                merged=merged.merge(calls,on=["chrom","start","end"],how="left",validate="one_to_one")
                merged["feature_id"]=[feature_id("insulation",c["REFERENCE_ID"],x,c["INSULATION_RESOLUTION"]) for x in merged[["chrom","start","end"]].itertuples(index=False,name=None)]
                model_input=directory/"insulation_model_scores.tsv";save(merged,model_input);outputs.append(model_input)
                if comparison["inference"]!="descriptive":run_r(c,"differential_insulation.R",[model_input,params,output])
                else:
                    result=merged[["feature_id","chrom","start","end"]].copy();result["case_minus_reference_log2_insulation"]=merged[case].mean(axis=1)-merged[control].mean(axis=1)
                    result["PValue"]=np.nan;result["FDR"]=np.nan;result["inference_status"]="descriptive";save(result,output)
            result=pd.read_csv(output,sep="\t")
            result["validation_status"]=status
            save(result,output)
            outputs.append(output)
    return outputs


def differential(c,scope):
    modules=set(c.get('_comparison_modules',[])) or {module for comparison in validate_metadata(c)['comparisons'] for module in comparison['modules'].split(';')};outputs=[]
    if "insulation" in modules:outputs.extend(_differential_single(dict(c,_comparison_modules=["insulation"]),scope))
    if "loops" in modules:
        for res in quantification_resolutions(c):
            products=_differential_single(dict(c,LOOP_RESOLUTION=res,_loop_measurement_resolution=res,_comparison_modules=["loops"]),scope)
            outputs.extend(products)
            if res==c["LOOP_RESOLUTION"]:
                for path in products:
                    alias=path.parent.parent/path.name
                    # Keep shared parameters when insulation already provided them.
                    if path.name=="parameters.json" and alias in outputs:continue
                    shutil.copyfile(path,alias);outputs.append(alias)
    return list(dict.fromkeys(outputs))


def depth_sensitivity(c,scope):
    root=Path(c["run_dir"]);samples=validate_metadata(c)["samples"]
    resolutions=resolution_set(c);fine=min(resolutions)
    require(all(res%fine==0 for res in resolutions),"Sensitivity grids must be integer coarsenings of the finest requested resolution")
    depths={};contexts={}
    for sample in samples:
        sid=sample["sample_id"];contexts[sid]=context(c,sid,fine)
        clr,bins,mask,labels,views,values=contexts[sid]
        cutoff=read_json(sample["provenance_path"])["filters"]["min_distance_bp"]
        total=0;positions=bins.start.to_numpy()
        for chunk in pixel_chunks(clr,c["PIXEL_CHUNKSIZE"]):
            i,j,y=[chunk[k].to_numpy() for k in ("bin1_id","bin2_id","count")]
            valid=mask[i]&mask[j]&(labels[i]>=0)&(labels[i]==labels[j])&(positions[j]-positions[i]>=cutoff)
            total+=int(y[valid].sum())
        require(total>0,"Cannot downsample an empty eligible cis domain");depths[sid]=total
    target=min(depths.values());rows=[];outputs=[]
    for seed_index in range(c["DOWNSAMPLING_REPEATS"]):
        seed_maps=[]
        for sample_index,sample in enumerate(samples):
            sid=sample["sample_id"];clr,bins,mask,labels,views,values=contexts[sid]
            common_bins=pd.read_csv(root/f"00_metadata/masks/{fine}.bins.tsv",sep="\t")
            ids=native_indices(clr,common_bins)
            mapping=np.full(len(bins),-1,dtype=np.int64);mapping[ids]=np.arange(len(common_bins))
            rng=np.random.default_rng(np.random.SeedSequence([c["RANDOM_SEED"],seed_index,sample_index]))
            probability=target/depths[sid];positions=bins.start.to_numpy()
            cutoff=read_json(sample["provenance_path"])["filters"]["min_distance_bp"]
            directory=root/f"01_qc/depth_sensitivity/seed_{seed_index}/{sid}";directory.mkdir(parents=True,exist_ok=True)
            fine_path=directory/f"{fine}.cool"
            def draw():
                for chunk in pixel_chunks(clr,c["PIXEL_CHUNKSIZE"]):
                    i,j=[chunk[k].to_numpy() for k in ("bin1_id","bin2_id")]
                    eligible=mask[i]&mask[j]&(labels[i]>=0)&(labels[i]==labels[j])&(positions[j]-positions[i]>=cutoff)
                    subset=chunk.loc[eligible].copy()
                    a=mapping[subset.bin1_id.to_numpy()];b=mapping[subset.bin2_id.to_numpy()]
                    subset["bin1_id"]=np.minimum(a,b);subset["bin2_id"]=np.maximum(a,b)
                    subset["count"]=rng.binomial(subset["count"].to_numpy(dtype=np.int64),probability)
                    yield subset.loc[subset["count"]>0].sort_values(["bin1_id","bin2_id"])
            scratch=Path(c["SCRATCH_ROOT"]);scratch.mkdir(parents=True,exist_ok=True)
            cooler.create_cooler(str(fine_path),common_bins,draw(),ordered=False,
                                 mergebuf=c["PIXEL_CHUNKSIZE"],temp_dir=str(scratch),
                                 dtypes={"count":"int64"},assembly=c["REFERENCE_ID"])
            realized=int(cooler.Cooler(str(fine_path)).info["sum"])
            rows.append({"sample_id":sid,"seed_index":seed_index,"base_seed":c["RANDOM_SEED"],"probability":probability,"native_eligible_depth":depths[sid],"target_expected_depth":target,"realized_depth":realized,"generator":"numpy.PCG64","numpy_version":np.__version__,"inference_status":"descriptive_technical_sensitivity"})
            maps=[fine_path]
            for res in resolutions:
                if res==fine:continue
                path=directory/f"{res}.cool"
                cooler.coarsen_cooler(str(fine_path),str(path),factor=res//fine,chunksize=c["PIXEL_CHUNKSIZE"],nproc=1)
                require(int(cooler.Cooler(str(path)).info["sum"])==realized,"Thinned coarsening violates count conservation");maps.append(path)
            for path in maps:
                raw=cooler.Cooler(str(path));res=raw.binsize
                archive=np.load(root/f"00_metadata/masks/{res}.npz");localmask=archive["mask"]
                localbins=raw.bins()[:][["chrom","start","end"]];localviews=pd.read_csv(root/f"00_metadata/masks/{res}.views.tsv",sep="\t")
                localvalues=np.ones(len(localbins))
                if c["DEPTH_SENSITIVITY_MODE"]=="balanced_features":
                    original=cooler.Cooler(matrix_uri(sample_matrix(c,sample,res),res))
                    unchanged=same_raw_counts(raw,original,c["PIXEL_CHUNKSIZE"])
                    if unchanged:localvalues=archive[sid].copy()
                    else:
                        balance_new_map(c,path,blacklist=np.flatnonzero(~localmask),cis_only=True)
                        localvalues=raw.bins()["weight"][:].to_numpy()
                    localmask=localmask&np.isfinite(localvalues)
                localvalues[~localmask]=np.nan
                expected=expected_table(raw,localbins,localmask,archive["labels"],localviews,localvalues,c)
                expected["normalization"]="raw_only" if c["DEPTH_SENSITIVITY_MODE"]=="raw_only" else "new_cis_only_derived_weights"
                if c["DEPTH_SENSITIVITY_MODE"]=="raw_only":expected[["balanced.sum","balanced.avg"]]=np.nan
                elif unchanged:expected["normalization"]="inherited_weights_for_unchanged_counts"
                outputs.append(save(expected,directory/f"{res}.expected.tsv"))
                if c["DEPTH_SENSITIVITY_MODE"]=="balanced_features" and res==c["INSULATION_RESOLUTION"] and c["RUN_INSULATION"]:
                    import cooltools
                    masked=MaskedCooler(str(path),localvalues)
                    table=cooltools.insulation(masked,[c["INSULATION_PRIMARY_WINDOW_BP"]],view_df=tool_views(masked,localviews),ignore_diags=c["INSULATION_IGNORE_DIAGS"],min_frac_valid_pixels=c["INSULATION_MIN_VALID_PIXEL_FRACTION"],nproc=1)
                    table["inference_status"]="descriptive_technical_sensitivity"
                    outputs.append(save(table,directory/"insulation.tsv"))
                if c["DEPTH_SENSITIVITY_MODE"]=="balanced_features" and res==c["LOOP_RESOLUTION"]:
                    feature_path=root/"04_quantification/discovered_loops/counts.tsv"
                    if feature_path.exists():
                        fixed=pd.read_csv(feature_path,sep="\t")[LOOP_COLUMNS]
                        counts,exposure,reasons=quantify(raw,fixed,localbins,localmask,archive["labels"],c)
                        summary=fixed.copy();summary["raw_count"]=counts;summary["relative_allocation"]=counts/exposure
                        exp=expected.set_index(["region1","dist"])["balanced.avg"]
                        lookup={(str(b.chrom),int(b.start)):i for i,b in enumerate(localbins.itertuples())}
                        oe=[]
                        for k,row in fixed.iterrows():
                            i=lookup.get((row.chrom1,row.start1));j=lookup.get((row.chrom2,row.start2))
                            if k in reasons or i is None or j is None:oe.append(np.nan);continue
                            view=localviews.iloc[archive["labels"][i]]["name"]
                            e=exp.get((view,j-i),np.nan)
                            oe.append(float(counts.iloc[k])*localvalues[i]*localvalues[j]/e if np.isfinite(e) and e>0 else np.nan)
                        summary["balanced_observed_expected"]=oe;summary["inference_status"]="descriptive_technical_sensitivity"
                        outputs.append(save(summary,directory/"fixed_loops.tsv"))
                if c["KEEP_DOWNSAMPLED_MAPS"]:outputs.append(path)
                seed_maps.append(path)
        if c["RUN_HICREP"]:
            config_path=root/"00_metadata/config.json"
            helper=str(Path(__file__).parents[2]/"scripts/scientific.py")
            command=([c["PYTHON_BOOTSTRAP"]] if c["ENVIRONMENT_MODE"]=="development" else [c["CONDA_EXE"],"run","--no-capture-output","-p",c["HICREP_ENV_PREFIX"],"python"])
            derived=root/f"01_qc/depth_sensitivity/seed_{seed_index}"
            subprocess.run([*command,helper,"depth_hicrep","--config",str(config_path),"--derived-root",str(derived)],check=True)
            outputs.extend([derived/"01_qc/hicrep/by_view.tsv",derived/"01_qc/hicrep/summary.tsv"])
        if not c["KEEP_DOWNSAMPLED_MAPS"]:
            for path in seed_maps:
                require(path.resolve().is_relative_to((root/"01_qc/depth_sensitivity").resolve()),"Thinned-map cleanup target escapes task root")
                path.unlink()
    output=root/"01_qc/depth_sensitivity/summary.tsv";outputs.append(save(pd.DataFrame(rows),output))
    return outputs


def cleanup_derived(c):
    """Discard only registered task-owned maps after all consumers succeeded."""
    root=Path(c["run_dir"]).resolve();registry_path=root/"02_pools/registry.json"
    if not c["RUN_POOL"] or c["KEEP_POOLED_MAPS"] or not registry_path.exists():return
    registry=read_json(registry_path)
    removed=set()
    for pool in registry["pools"]:
        if not pool.get("task_owned"):continue
        for uri in pool["uris"].values():
            path=Path(uri.split("::")[0]).resolve()
            require(path.is_relative_to(root/"02_pools"),"Cleanup target escapes task-owned pools directory")
            if path.exists():path.unlink();removed.add(str(path.relative_to(root)))
    if removed:
        from .state import record_path
        record=read_json(record_path(c,"pool","run"))
        for output in record["outputs"]:
            if output["path"] in removed:output["discarded_intermediate"]=True
        atomic_json(record_path(c,"pool","run"),record)


IMPLEMENTATIONS={"preflight":lambda c,s:inspect_inputs(c,writing=True),"register":register,"qc":qc,
                 "pool":pools,"compartments":compartments,"insulation":insulation,"tads":lambda c,s:__import__("microc2compare.tads",fromlist=["run"]).run(c,s),"loops":loops,
                 "consensus":consensus,"quantify":quantify_features,"aggregates":aggregates,
                 "chromosight_discovery":lambda c,s:__import__("microc2compare.fine_loops",fromlist=["chromosight_discovery"]).chromosight_discovery(c,s),
                 "stripes":lambda c,s:__import__("microc2compare.stripes",fromlist=["stripes"]).stripes(c,s),
                 "stripe_consensus":lambda c,s:__import__("microc2compare.stripes",fromlist=["stripes_summary"]).stripes_summary(c),
                 "browser_export":lambda c,s:__import__("microc2compare.browser_export",fromlist=["export_browser"]).export_browser(c,s),
                 "chromosight":lambda c,s:__import__("microc2compare.fine_loops",fromlist=["chromosight_scores"]).chromosight_scores(c,s),
                 "annotation":lambda c,s:__import__("microc2compare.annotation",fromlist=["run"]).run(c,s),
                 "differential":differential,"depth_sensitivity":depth_sensitivity,"hicrep":hicrep}
