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
                       consensus_loops, feature_id, quantify)
from .matrix import (MaskedCooler, bounded_raw, domain, expected_table, matrix_uri,
                     pixel_chunks, raw_summary, reference, resolution_set, unordered_total, weights)
from .state import inventory, source_digest, scopes as stage_scopes
from .analyses import aggregate_modules, comparison_modules


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
    require(registry["reference_id"] == c["REFERENCE_ID"] and registry["resolution"] == c["LOOP_RESOLUTION"], "External feature reference/resolution mismatch")
    path = resolve(registry["path"], Path(c["INDEPENDENT_FEATURES_REGISTRY"]).parent)
    require(file_hash(path) == registry["sha256"], "Frozen external feature checksum mismatch")
    require(not set(registry["discovery_samples"]) & {s["sample_id"] for s in samples}, "External loop discovery reuses analysis samples")
    require(registry["independent_of_analysis_samples"] is True, "Formal loop features must be independently discovered")
    require(registry.get("source") and registry.get("independence_evidence"), "External candidates require inspectable source/independence evidence")
    return canonical_loops(pd.read_csv(path, sep="\t"), c["REFERENCE_ID"], c["LOOP_RESOLUTION"]), registry


def environment_record(c, kind):
    if c["ENVIRONMENT_MODE"] == "development":
        if kind == "core": return inventory()
        if kind == "hicrep": return {"hicrep": inventory()["packages"]["hicrep"]}
        command = [c["R_SCRIPT"], "--vanilla", "-e", 'cat(R.version.string,"\\n"); for(p in c("edgeR","limma","statmod")) cat(p,as.character(packageVersion(p)),"\\n")']
    else:
        prefix = c[{"core": "CORE_ENV_PREFIX", "hicrep": "HICREP_ENV_PREFIX", "diff_r": "DIFF_R_ENV_PREFIX"}[kind]]
        require((Path(prefix)/"conda-meta").is_dir(), f"Missing {kind} prefix {prefix}; run explicit bootstrap")
        packages = subprocess.run([c["CONDA_EXE"], "list", "-p", prefix, "--json"], check=True, text=True, capture_output=True)
        if c["ENVIRONMENT_MODE"] == "locked":
            lock = c[{"core":"CORE_LOCK","hicrep":"HICREP_LOCK","diff_r":"DIFF_R_LOCK"}[kind]]
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
    keys = ("LOOP_RESOLUTION", "LOOP_MIN_DISTANCE_BP", "LOOP_MAX_DISTANCE_BP", "LOOP_MIN_TOTAL_COUNT", "LOOP_MIN_MEAN_COUNT_AT_MIN_EXPOSURE", "LOOP_MIN_FEATURES_FORMAL", "INSULATION_RESOLUTION", "INSULATION_PRIMARY_WINDOW_BP", "INSULATION_MIN_VALID_PIXEL_FRACTION", "FDR_THRESHOLD")
    if design_domain is None:design_domain=calibration_design_domain(c,comparison,validate_metadata(c)["samples"])
    return digest({"source": source_digest("differential"), "environment": environments,"design_domain":design_domain,
                   "module": module, "design": formula_terms(comparison["design"]),
                   "recipe": {k:c[k] for k in keys}})


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
    for s in metadata["samples"]:
        prov = read_json(s["provenance_path"])
        jsonschema.validate(prov, read_json(root/"schemas/provenance.schema.json"))
        own = [r for r in metadata["membership"] if r["sample_id"] == s["sample_id"]]
        require(set(prov["library_ids"]) == {r["library_id"] for r in own} and set(prov["source_ids"]) == {r["source_id"] for r in own}, f"{s['sample_id']}: upstream membership does not match biological manifest")
        require(c["LOOP_MIN_DISTANCE_BP"] >= prov["filters"]["min_distance_bp"], "Loop distance starts below upstream contact cutoff")
    for res in resolution_set(c):
        domain(c, metadata["samples"], res)
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
    feature_set, feature_registry = external_features(c, metadata["samples"]) if "loops" in comparison_modules(c) or "apa" in aggregate_modules(c) or c["RUN_LOOPS"] else (pd.DataFrame(columns=LOOP_COLUMNS), None)
    environments = {"core": environment_record(c, "core")}
    needed_tools = c["RUN_COMPARTMENTS"] or c["RUN_INSULATION"] or c["RUN_LOOPS"]
    if needed_tools:
        import cooltools
        require(hasattr(cooltools, "eigs_cis") and hasattr(cooltools, "insulation") and hasattr(cooltools, "dots"), "Core cooltools API incompatibility")
    if c["RUN_HICREP"]:
        environments["hicrep"] = environment_record(c, "hicrep")
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
            check_calibration(c, comparison, environments)
    if writing:
        out = Path(c["run_dir"])/"00_metadata"
        prior = read_json(out/"inputs.json")
        for key in ("chromsizes", "views", "blacklist", "phasing_track"):
            if key == "phasing_track" and not c["RUN_COMPARTMENTS"]: continue
            if ref.get(key): prior[ref[key]] = fingerprint(ref[key], prior.get(ref[key]))
        if feature_registry:
            path = resolve(feature_registry["path"], Path(c["INDEPENDENT_FEATURES_REGISTRY"]).parent)
            prior[path] = fingerprint(path, prior.get(path))
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
    return outputs


def context(c, scope, resolution):
    base = Path(c["run_dir"])/f"00_metadata/masks/{resolution}"
    archive = np.load(str(base)+".npz")
    bins = pd.read_csv(str(base)+".bins.tsv", sep="\t")
    views = pd.read_csv(str(base)+".views.tsv", sep="\t")
    if scope.startswith("pool__"):
        pool = next(p for p in read_json(Path(c["run_dir"])/"02_pools/registry.json")["pools"] if p["pool_id"] == scope[6:])
        uri = pool["uris"][str(resolution)]
        raw = cooler.Cooler(uri)
        pool_config=dict(c,WEIGHT_COLUMN="weight",WEIGHT_CONVENTION="multiplicative") if pool.get("task_owned") else c
        values = weights(raw, pool_config, read_json(pool["provenance_path"]), resolution)
        values[~archive["mask"]] = np.nan
    else:
        sample = next(s for s in validate_metadata(c)["samples"] if s["sample_id"] == scope)
        uri = matrix_uri(sample["matrix_path"], resolution)
        values = archive[scope]
    mask = archive["mask"] & np.isfinite(values)
    return MaskedCooler(uri, values), bins, mask, archive["labels"], views, values


def qc(c, scope):
    outputs = []
    if not c["RUN_QC"]:
        path=Path(c["run_dir"])/f"01_qc/hicrep/inputs/{scope}.json"
        atomic_json(path,{"sample_id":scope,"status":"raw_only_HiCRep_input_validated","resolution":c["HICREP_RESOLUTION"],"normalization":"raw_counts_no_balancing"})
        return [path]
    all_res={c["PS_RESOLUTION"]}
    if c["RUN_LOOPS"]:all_res.add(c["LOOP_RESOLUTION"])
    if c["RUN_LOOPS"] and c["LOOP_SENSITIVITY_RESOLUTIONS"]:all_res.add(int(c["LOOP_SENSITIVITY_RESOLUTIONS"]))
    all_res=sorted(all_res)
    for res in all_res:
        clr, bins, mask, labels, views, values = context(c, scope, res)
        summary = raw_summary(clr, bins, mask, labels, c)
        path = Path(c["run_dir"])/f"01_qc/contacts/{scope}.{res}.json"
        atomic_json(path, dict(summary, sample_id=scope, resolution=res, normalization="inherited_upstream_weights"))
        expected = expected_table(clr, bins, mask, labels, views, values, c)
        outputs += [path, save(expected, Path(c["run_dir"])/f"01_qc/ps/{scope}.{res}.tsv")]
    return outputs


def phasing_track(c, bins):
    ref, _, _ = reference(c)
    require(ref.get("phasing_track"), "Compartments require a registered fixed phasing track")
    track = pd.read_csv(ref["phasing_track"], sep="\t", header=None, names=["chrom", "start", "end", "value"])
    # Exact grid is required; no hidden interpolation/lift-over.
    merged = bins.merge(track, on=["chrom", "start", "end"], how="left", validate="one_to_one")
    require(merged.value.notna().any(), "Phasing track has no matching bins")
    return merged


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
        eigenvalues, vectors = cooltools.eigs_cis(clr, view_df=views, n_eigs=3, clr_weight_name="weight", ignore_diags=c["INSULATION_IGNORE_DIAGS"],map=pool.map)
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
    table = cooltools.insulation(clr, windows, view_df=views, clr_weight_name="weight", ignore_diags=c["INSULATION_IGNORE_DIAGS"], min_frac_valid_pixels=c["INSULATION_MIN_VALID_PIXEL_FRACTION"], threshold=c["INSULATION_BOUNDARY_THRESHOLD"], append_raw_scores=True, chunksize=c["PIXEL_CHUNKSIZE"], nproc=c["THREADS_INSULATION"])
    for window in windows:
        width=window//c["INSULATION_RESOLUTION"]
        theoretical=sum(a+b>=c["INSULATION_IGNORE_DIAGS"] for a in range(width) for b in range(width))
        fraction=table[f"n_valid_pixels_{window}"]/theoretical if theoretical else np.zeros(len(table))
        table[f"valid_pixel_fraction_{window}"]=fraction
        insufficient=fraction<c["INSULATION_MIN_VALID_PIXEL_FRACTION"]
        table.loc[insufficient,[f"log2_insulation_score_{window}",f"boundary_strength_{window}"]]=np.nan
        table.loc[insufficient,f"is_boundary_{window}"]=False
    if "region" not in table:
        table["region"] = [views.iloc[labels[i]]["name"] if labels[i] >= 0 else pd.NA for i in table.index]
    table["inference_status"] = "per_sample_measurement"
    outputs = [save(table, Path(c["run_dir"])/f"03_features/insulation/{scope}.tsv")]
    for window in windows:
        outputs.append(track(table,f"log2_insulation_score_{window}",Path(c["run_dir"])/f"03_features/insulation/{scope}.{window}.bedgraph"))
        boundary = table.loc[table[f"is_boundary_{window}"].fillna(False), ["chrom", "start", "end", "region", f"boundary_strength_{window}"]]
        outputs += [save(boundary, Path(c["run_dir"])/f"03_features/boundaries/{scope}.{window}.tsv")]
    return outputs


def loops(c, scope):
    import cooltools
    outputs = []
    for res in [c["LOOP_RESOLUTION"]] + ([int(c["LOOP_SENSITIVITY_RESOLUTIONS"])] if c["LOOP_SENSITIVITY_RESOLUTIONS"] else []):
        clr, bins, mask, labels, views, values = context(c, scope, res)
        expected_path = Path(c["run_dir"])/f"01_qc/ps/{scope}.{res}.tsv"
        expected = pd.read_csv(expected_path, sep="\t") if expected_path.exists() else expected_table(clr, bins, mask, labels, views, values, c)
        calls = cooltools.dots(clr, expected, expected_value_col="balanced.avg", clr_weight_name="weight", view_df=views, max_loci_separation=c["LOOP_MAX_DISTANCE_BP"], nproc=c["THREADS_LOOPS"])
        if len(calls):
            calls = calls.loc[(calls.chrom1 == calls.chrom2) & (calls.start2-calls.start1 >= c["LOOP_MIN_DISTANCE_BP"]) & (calls.start2-calls.start1 <= c["LOOP_MAX_DISTANCE_BP"])]
            calls = canonical_loops(calls, c["REFERENCE_ID"], res)
        else:
            calls = pd.DataFrame(columns=LOOP_COLUMNS)
        calls["caller"] = "cooltools_dots"; calls["resolution"] = res
        calls["inference_status"] = "discovery_enrichment_only"
        outputs += [save(calls, Path(c["run_dir"])/f"03_features/loops/{scope}.{res}.bedpe.tsv")]
    return outputs


def consensus(c, scope):
    root = Path(c["run_dir"])
    outputs = []
    frames = {sid:pd.read_csv(root/f"03_features/loops/{sid}.{c['LOOP_RESOLUTION']}.bedpe.tsv",sep="\t") for sid in stage_scopes(c,"loops")} if c["RUN_LOOPS"] else {}
    loops_frame = consensus_loops(frames, c["REFERENCE_ID"], c["LOOP_RESOLUTION"], c["LOOP_CONSENSUS_TOLERANCE_BINS"])
    outputs += [save(loops_frame, root/"03_features/loops/consensus.tsv")]
    domains = []
    if c["RUN_INSULATION"]:
        for window in map(int, c["INSULATION_WINDOWS_BP"].split(",")):
            boundaries = {sid:pd.read_csv(root/f"03_features/boundaries/{sid}.{window}.tsv",sep="\t") for sid in stage_scopes(c,"insulation")}
            merged = boundary_consensus(boundaries, c["REFERENCE_ID"], c["INSULATION_RESOLUTION"], window)
            outputs += [save(merged, root/f"03_features/boundaries/consensus.{window}.tsv")]
            archive = np.load(root/f"00_metadata/masks/{c['INSULATION_RESOLUTION']}.npz")
            bins = pd.read_csv(root/f"00_metadata/masks/{c['INSULATION_RESOLUTION']}.bins.tsv", sep="\t")
            for (chrom, region), group in merged.groupby(["chrom", "region"]):
                positions = sorted(group.start)
                for start, end in zip(positions, positions[1:]):
                    inside = (bins.chrom.astype(str) == chrom) & (bins.start >= start) & (bins.end <= end)
                    if not inside.any() or not archive["mask"][inside].all(): continue
                    domains.append({"chrom":chrom, "start":start, "end":end, "region":region, "window_bp":window,
                                    "feature_id":feature_id("domain",c["REFERENCE_ID"],(chrom,start,end,window),c["INSULATION_RESOLUTION"]), "inference_status":"descriptive"})
    outputs += [save(pd.DataFrame(domains, columns=["chrom","start","end","region","window_bp","feature_id","inference_status"]), root/"03_features/boundary_derived_domains/consensus.tsv")]
    return outputs


def quantify_features(c, scope):
    root = Path(c["run_dir"])
    samples = validate_metadata(c)["samples"]
    discovered_path = root/"03_features/loops/consensus.tsv"
    discovered = pd.read_csv(discovered_path,sep="\t") if c["RUN_LOOPS"] and discovered_path.exists() else pd.DataFrame(columns=LOOP_COLUMNS)
    independent, registry = external_features(c,samples)
    outputs=[]
    for name, universe in (("independent_loops",independent),("discovered_loops",discovered)):
        universe = universe.reset_index(drop=True)
        matrix=universe[LOOP_COLUMNS].copy()
        exposure_rows, exclusions=[],[]
        for sample in samples:
            clr,bins,mask,labels,views,values=context(c,sample["sample_id"],c["LOOP_RESOLUTION"])
            counts,exposure,reasons=quantify(clr,universe,bins,mask,labels,c)
            matrix[sample["sample_id"]]=counts
            exposure_rows.append({"sample_id":sample["sample_id"],"eligible_cis_exposure":exposure})
            exclusions.extend({"sample_id":sample["sample_id"],"feature_id":universe.loc[i,"feature_id"],"reason":reason} for i,reason in reasons.items())
        directory=root/"04_quantification"/name
        outputs += [save(matrix,directory/"counts.tsv"),save(pd.DataFrame(exposure_rows),directory/"exposures.tsv"),save(pd.DataFrame(exclusions,columns=["sample_id","feature_id","reason"]),directory/"exclusions.tsv")]
        provenance=directory/"features.json"
        atomic_json(provenance,{"inference_status":"independent_candidates" if name=="independent_loops" and registry else "descriptive_same_data", "registry":registry if name=="independent_loops" else None, "feature_hash":digest(universe[LOOP_COLUMNS].to_dict("records"))})
        outputs.append(provenance)
    return outputs


def oe_window(clr, region1, region2, values, labels, views, expected, c):
    lo,hi=clr.extent(region1); lo2,hi2=clr.extent(region2)
    raw=bounded_raw(clr,region1,region2,c).astype(float)
    result=raw*values[lo:hi,None]*values[None,lo2:hi2]
    i=np.arange(lo,hi)[:,None]; j=np.arange(lo2,hi2)[None,:]
    distance=np.abs(j-i)
    denominator=np.full(result.shape,np.nan)
    for vid in np.unique(labels[lo:hi]):
        if vid<0:continue
        name=views.iloc[vid]["name"]
        ex=expected.loc[(expected.region1==name)&(expected.region2==name)].set_index("dist")["balanced.avg"]
        same=(labels[lo:hi,None]==vid)&(labels[None,lo2:hi2]==vid)
        lookup=np.full(max(int(distance.max())+1,1),np.nan)
        indices=ex.index.to_numpy(dtype=int)
        valid=indices<len(lookup)
        lookup[indices[valid]]=ex.to_numpy()[valid]
        denominator[same]=lookup[distance[same]]
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
        ids=["pool__"+p["pool_id"] for p in read_json(root/"02_pools/registry.json")["pools"] if p["pool_scope"]=="all_samples"] if c["SADDLE_PARTITION_SOURCE"]=="all_sample_pool" else [s["sample_id"] for s in samples]
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


def aggregates(c,scope):
    modules=aggregate_modules(c)
    root=Path(c["run_dir"])
    samples=validate_metadata(c)["samples"]
    loops_path=root/"03_features/loops/consensus.tsv"
    discovered=pd.read_csv(loops_path,sep="\t") if "apa" in modules and c["RUN_LOOPS"] and loops_path.exists() else pd.DataFrame(columns=LOOP_COLUMNS)
    independent,_=external_features(c,samples) if "apa" in modules else (pd.DataFrame(columns=LOOP_COLUMNS), None)
    universe=independent if len(independent) else discovered
    domains_path=root/"03_features/boundary_derived_domains/consensus.tsv"
    domains=pd.read_csv(domains_path,sep="\t") if "ada" in modules and c["RUN_INSULATION"] and domains_path.exists() else pd.DataFrame()
    if len(domains):domains=domains.loc[domains.window_bp==c["INSULATION_PRIMARY_WINDOW_BP"]]
    outputs=[]; summaries=[]; domain_rows=[]
    for sample in samples if modules & {"apa", "ada"} else []:
        sid=sample["sample_id"];res=c["LOOP_RESOLUTION"]
        clr,bins,mask,labels,views,values=context(c,sid,res)
        expected=expected_table(clr,bins,mask,labels,views,values,c)
        if "apa" in modules:
            flank=100000//res; shape=2*flank+1
            total=np.zeros((shape,shape));coverage=np.zeros((shape,shape),dtype=int)
            for feature in universe.itertuples():
                region1=(feature.chrom1,feature.start1-100000,feature.end1+100000)
                region2=(feature.chrom2,feature.start2-100000,feature.end2+100000)
                if region1[1]<0 or region2[1]<0 or region1[2]>clr.chromsizes[feature.chrom1] or region2[2]>clr.chromsizes[feature.chrom2]:continue
                window=oe_window(clr,region1,region2,values,labels,views,expected,c)
                if window.shape!=total.shape:continue
                good=np.isfinite(window)
                total[good]+=window[good];coverage[good]+=1
            average=np.divide(total,coverage,out=np.full(total.shape,np.nan),where=coverage>0)
            path=root/f"05_aggregates/apa/{sid}.npz";path.parent.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(path,mean=average,sum=total,valid_feature_count=coverage)
            outputs.append(path)
            summaries.append({"sample_id":sid,"central_oe":average[flank,flank],"central_valid_features":coverage[flank,flank],"central_3x3_oe":float(np.nanmean(average[flank-1:flank+2,flank-1:flank+2])) if np.isfinite(average).any() else np.nan,"inference_status":"descriptive"})
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
            np.savez_compressed(path,mean=np.divide(ada_total,ada_count,out=np.full((50,50),np.nan),where=ada_count>0),sum=ada_total,valid_domain_count=ada_count)
            outputs.append(path)
    if "apa" in modules:
        outputs.append(save(pd.DataFrame(summaries),root/"05_aggregates/apa/summary.tsv"))
    if "ada" in modules:
        outputs.append(save(pd.DataFrame(domain_rows,columns=["feature_id","sample_id","mean_intra_domain_oe","oe_sum","eligible_pixel_count","inference_status"]),root/"04_quantification/domain_scores/strength.tsv"))
    if "saddle" in modules:
        partition,labels=saddle_partition(c,samples)
        outputs.append(save(partition,root/"05_aggregates/saddle/partition.tsv"))
        info=root/"05_aggregates/saddle/partition.json"
        atomic_json(info,{"source":c["SADDLE_PARTITION_SOURCE"],"groups":50,"partition_hash":digest(partition[["chrom","start","end","partition"]].to_dict("records")),"views":[{"view_index":int(vid),"usable_bins":int(((labels==vid)&(partition.partition.to_numpy()>=0)).sum()),"status":"partitioned" if ((labels==vid)&(partition.partition.to_numpy()>=0)).sum()>=50 else "fewer_than_50_usable_bins"} for vid in np.unique(labels[labels>=0])],"inference_status":"descriptive"})
        outputs.append(info)
        saddle_rows=[]
        for sample in samples:
            sid=sample["sample_id"];clr,bins,mask,labels,views,values=context(c,sid,c["COMPARTMENT_RESOLUTION"])
            ex=expected_table(clr,bins,mask,labels,views,values,c)
            assignments=partition.partition.to_numpy();sums=np.zeros((50,50));counts=np.zeros((50,50),dtype=np.int64)
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
            np.savez_compressed(path,mean=average,sum=sums,valid_pixel_count=counts);outputs.append(path)
            def corner(x,y):
                n=counts[x,y].sum();return sums[x,y].sum()/n if n else np.nan
            aa=corner(slice(40,50),slice(40,50));bb=corner(slice(0,10),slice(0,10));ab=corner(slice(0,10),slice(40,50))
            saddle_rows.append({"sample_id":sid,"AA":aa,"BB":bb,"AB":ab,"segregation_ratio":(aa+bb)/(2*ab) if np.isfinite(ab) and ab>0 else np.nan,"inference_status":"descriptive"})
        outputs.append(save(pd.DataFrame(saddle_rows),root/"05_aggregates/saddle/summary.tsv"))
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
        source=cooler.Cooler(matrix_uri(sample["matrix_path"],clr.binsize))
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
                actual=raw.bins()[:][["chrom","start","end"]].copy();actual["chrom"]=actual.chrom.astype(str)
                require(np.array_equal(actual.to_numpy(),bins.to_numpy()),"Pool axes mismatch")
                weights(raw,c,prov,res)
                total=sum(int(p["count"].sum()) for p in pixel_chunks(raw,c["PIXEL_CHUNKSIZE"]))
                wanted=sum(unordered_total(cooler.Cooler(matrix_uri(s["matrix_path"],res)),c["PIXEL_CHUNKSIZE"]) for s in samples if s["sample_id"] in members)
                require(total==wanted,"Imported pool raw total differs from member sum")
                validate_pool_sum(raw,[cooler.Cooler(matrix_uri(s["matrix_path"],res)) for s in samples if s["sample_id"] in members],c["PIXEL_CHUNKSIZE"])
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
                registered.append({"pool_id":pool_id,"pool_scope":kind,"condition":condition,"member_sample_ids":sample["sample_id"],"provenance_path":sample["provenance_path"],"task_owned":False,"uris":{str(res):matrix_uri(sample["matrix_path"],res) for res in resolution_set(c)}})
                continue
            require(c["BALANCE_NEW_DERIVED_MAPS"],"Creating pools for balanced discovery requires explicit new-map balancing")
            row={"pool_id":pool_id,"pool_scope":kind,"condition":condition,"member_sample_ids":";".join(s["sample_id"] for s in members),"task_owned":True,"uris":{}}
            input_records=read_json(root/"00_metadata/inputs.json")
            source_prov=read_json(members[0]["provenance_path"])
            prov={"reference_id":c["REFERENCE_ID"],"count_units":"unique_contacts","filters":source_prov["filters"],"balancing":{},"member_sample_ids":[s["sample_id"] for s in members],"member_hashes":{s["sample_id"]:input_records[s["matrix_path"].split("::")[0]]["sha256"] for s in members},"dedup_scope":"upstream_biological_members_already_deduplicated","library_ids":[lib for s in members for lib in read_json(s["provenance_path"])["library_ids"]],"source_ids":[source for s in members for source in read_json(s["provenance_path"])["source_ids"]],"summation":"raw_integer_sum"}
            for res in resolution_set(c):
                path=root/f"02_pools/{pool_id}/{res}.cool";path.parent.mkdir(parents=True,exist_ok=True)
                uris=[matrix_uri(s["matrix_path"],res) for s in members]
                require(all(cooler.Cooler(u).info.get("storage-mode")=="symmetric-upper" for u in uris),
                        "Explicit new pooling requires symmetric-upper members; request upstream pooling for square storage")
                pending=path.with_name("."+path.name+".pending.cool")
                cooler.merge_coolers(str(pending),uris,mergebuf=c["PIXEL_CHUNKSIZE"],columns=["count"],dtypes={"count":"int64"})
                require(cooler.Cooler(str(pending)).info["sum"]==sum(cooler.Cooler(u).info["sum"] for u in uris),"New pool violates raw sum conservation")
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


def differential(c,scope):
    root=Path(c["run_dir"]);metadata=validate_metadata(c);outputs=[]
    for comparison in metadata["comparisons"]:
        if c.get("_comparison_modules"):
            comparison = dict(comparison, modules=";".join(m for m in comparison["modules"].split(";") if m in c["_comparison_modules"]))
            if not comparison["modules"]: continue
        design=design_matrix(metadata["samples"],comparison,c)
        ids=design["sample_ids"]
        case=[s["sample_id"] for s in metadata["samples"] if s["condition"]==comparison["case_condition"]]
        control=[s["sample_id"] for s in metadata["samples"] if s["condition"]==comparison["reference_condition"]]
        directory=root/"06_comparisons"/comparison["comparison_id"];directory.mkdir(parents=True,exist_ok=True)
        status=check_calibration(c,comparison,read_json(root/"00_metadata/environments.json"))
        parameters={"design":design,"validation_status":status,"inference":comparison["inference"],"fdr":c["FDR_THRESHOLD"],"min_total":c["LOOP_MIN_TOTAL_COUNT"],"min_mean_at_min_exposure":c["LOOP_MIN_MEAN_COUNT_AT_MIN_EXPOSURE"],"min_features":c["LOOP_MIN_FEATURES_FORMAL"]}
        params=directory/"parameters.json";atomic_json(params,parameters);outputs.append(params)
        for module in comparison["modules"].split(";"):
            output=directory/(module+".tsv")
            if module=="loops":
                kind="independent_loops" if comparison["inference"]!="descriptive" else "discovered_loops"
                counts=pd.read_csv(root/f"04_quantification/{kind}/counts.tsv",sep="\t")
                if comparison["inference"]=="descriptive" and not len(counts):counts=pd.read_csv(root/"04_quantification/independent_loops/counts.tsv",sep="\t")
                exposures=pd.read_csv(root/f"04_quantification/{kind}/exposures.tsv",sep="\t").set_index("sample_id").eligible_cis_exposure
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
                    subset["count"]=rng.binomial(subset["count"].to_numpy(dtype=np.int64),probability)
                    yield subset.loc[subset["count"]>0]
            cooler.create_cooler(str(fine_path),bins,draw(),ordered=True,dtypes={"count":"int64"},assembly=c["REFERENCE_ID"])
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
                    original=cooler.Cooler(matrix_uri(sample["matrix_path"],res))
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
                    table=cooltools.insulation(masked,[c["INSULATION_PRIMARY_WINDOW_BP"]],view_df=localviews,ignore_diags=c["INSULATION_IGNORE_DIAGS"],min_frac_valid_pixels=c["INSULATION_MIN_VALID_PIXEL_FRACTION"],nproc=1)
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
                 "differential":differential,"depth_sensitivity":depth_sensitivity,"hicrep":hicrep}
