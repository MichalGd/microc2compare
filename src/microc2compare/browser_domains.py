"""Export existing TAD and ADA domain results without opening contact matrices.

This helper is used by the unified UCSC exporter. It deliberately distinguishes
the ADA attempted domain universe from its unrecorded, pixel-specific pileup
membership. All coordinates and feature identifiers originate in saved tables.
"""
from __future__ import annotations

import hashlib
import math
import re
import statistics


COORDS = ("chrom", "start", "end")
TAD_FIELDS = ("feature_id", *COORDS, "window_bp", "mean_log2_insulation",
              "left_boundary_strength", "right_boundary_strength", "support",
              "n_samples", "support_fraction", "sample_ids")
ADA_FIELDS = ("feature_id", *COORDS, "window_bp")
STRENGTH_FIELDS = ("feature_id", "sample_id", "mean_intra_domain_oe", "oe_sum",
                   "eligible_pixel_count", "inference_status")


def _slug(value):
    text = str(value)
    clean = re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("._") or "unnamed"
    return clean if clean == text else clean + "_" + hashlib.sha256(text.encode()).hexdigest()[:8]


def _number(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _integer(value, what):
    result = _number(value)
    if result is None or int(result) != result:
        raise ValueError(f"Noninteger {what}: {value!r}")
    return int(result)


def _coordinate(row):
    return str(row["chrom"]), _integer(row["start"], "start"), _integer(row["end"], "end")


def _required(ctx, source, fields):
    if not ctx.exists(source):
        raise ValueError(f"Partial domain module: missing {source}")
    return ctx.table(source, required=fields)


def _unique_ids(rows, source):
    result = {}
    for row in rows:
        fid = str(row["feature_id"])
        if not fid or fid in result:
            raise ValueError(f"Empty or duplicate domain feature_id in {source}: {fid!r}")
        result[fid] = row
    return result


def _bed9(row, color, identifier=None):
    chrom, start, end = _coordinate(row)
    return chrom, start, end, identifier or row["feature_id"], 0, ".", start, end, color


def _nonoverlapping(rows):
    previous = {}
    for row in sorted(rows, key=lambda x: _coordinate(x)):
        chrom, start, end = _coordinate(row)
        if start < previous.get(chrom, -1):
            return False
        previous[chrom] = end
    return True


def _tad_set(ctx, source, output, window, scope, label, supported=True):
    rows = _required(ctx, source, TAD_FIELDS)
    _unique_ids(rows, source)
    if any(_integer(row["window_bp"], "window_bp") != window for row in rows):
        raise ValueError(f"TAD window mismatch in {source}")
    ctx.copy_source(source, output + ".source.tsv")
    kind = "supported TADs" if supported else "unsupported or conflicting TAD candidates"
    ctx.bed(rel=output + ".bed", rows=[_bed9(row, "0,114,178" if supported else "128,128,128") for row in rows],
            description=f"{label}: {kind}; {window} bp insulation window; exact saved coordinates",
            analysis="tads", scope=scope, label=label, source=source,
            resolution=ctx.config.get("INSULATION_RESOLUTION", ""), window=window)
    finite = [row for row in rows if _number(row["mean_log2_insulation"]) is not None]
    if supported and finite:
        if not _nonoverlapping(finite):
            raise ValueError(f"Supported TAD intervals overlap in {source}; cannot export as bedGraph")
        ctx.bedgraph(rel=output + ".mean_log2_insulation.bedgraph",
                     rows=[(*_coordinate(row), float(row["mean_log2_insulation"])) for row in finite],
                     description=f"{label}: saved mean log2 insulation over each TAD; window {window} bp",
                     analysis="tads", scope=scope, label=label, source=source,
                     resolution=ctx.config.get("INSULATION_RESOLUTION", ""), window=window)
    return rows


def _export_tads(ctx):
    params_source = "03_features/tads/parameters.json"
    windows = [int(x) for x in str(ctx.config.get("INSULATION_WINDOWS_BP", "")).split(",") if x]
    if not ctx.exists(params_source):
        any_expected = any(ctx.exists(f"03_features/tads/samples/{sample['sample_id']}.{window}.tsv")
                           for sample in ctx.samples for window in windows)
        if any_expected:
            raise ValueError("TAD sample tables exist but 03_features/tads/parameters.json is missing")
        ctx.add_missing("TAD tracks", "The TAD result module is not available in this run copy.")
        return
    params = ctx.json(params_source)
    condition_ids = params.get("conditions", {})
    if set(condition_ids) != set(ctx.conditions):
        raise ValueError("TAD recorded conditions differ from the sample metadata")
    recorded = params.get("parameters", {})
    if not windows:
        raise ValueError("No recorded insulation windows for existing TAD module")
    for key in ("INSULATION_RESOLUTION", "INSULATION_WINDOWS_BP", "INSULATION_BOUNDARY_THRESHOLD"):
        if key in recorded and key in ctx.config and str(recorded[key]) != str(ctx.config[key]):
            raise ValueError(f"TAD parameter mismatch: {key}")
    ctx.copy_source(params_source, "tads/parameters.json")
    for window in windows:
        for sample in ctx.samples:
            sid = sample["sample_id"]
            source = f"03_features/tads/samples/{sid}.{window}.tsv"
            output = f"tads/individual/{_slug(sid)}/window_{window}"
            _tad_set(ctx, source, output, window, "individual", sid)
            exclusion_source = f"03_features/tads/samples/{sid}.{window}.exclusions.tsv"
            exclusions = _required(ctx, exclusion_source, (*COORDS, "window_bp", "reason"))
            ctx.copy_source(exclusion_source, output + ".exclusions.source.tsv")
            if any(_integer(row["window_bp"], "window_bp") != window for row in exclusions):
                raise ValueError(f"TAD exclusion window mismatch: {exclusion_source}")
            # Exclusions have no upstream IDs. Browser names explicitly say export-only.
            rows = []
            for row in exclusions:
                chrom, start, end = _coordinate(row)
                if start >= end:
                    raise ValueError(f"Invalid excluded TAD interval in {exclusion_source}")
                name = "excluded:" + row["reason"] + ":" + hashlib.sha256(
                    repr((sid, window, chrom, start, end, row["reason"])).encode()).hexdigest()[:12]
                rows.append(_bed9(row, "170,170,170", name))
            ctx.bed(rel=output + ".excluded_intervals.bed", rows=rows,
                    description=f"{sid}: excluded TAD intervals and recorded reasons; browser-only names",
                    analysis="tads", scope="individual", label=sid, source=exclusion_source,
                    resolution=ctx.config.get("INSULATION_RESOLUTION", ""), window=window)
        for condition, members in ctx.conditions.items():
            token = condition_ids[condition]
            source = f"03_features/tads/conditions/{token}.{window}.tsv"
            output = f"tads/conditions/{_slug(condition)}/window_{window}"
            supported = _tad_set(ctx, source, output, window, "condition", condition)
            unsupported_source = f"03_features/tads/conditions/{token}.{window}.unsupported.tsv"
            unsupported = _tad_set(ctx, unsupported_source, output + ".unsupported", window,
                                   "condition", condition, supported=False)
            for row in supported + unsupported:
                voters = set(filter(None, row["sample_ids"].split(";")))
                support = _integer(row["support"], "support")
                n = _integer(row["n_samples"], "n_samples")
                fraction = _number(row["support_fraction"])
                if not voters <= set(members) or len(voters) != support or n != len(members):
                    raise ValueError(f"Invalid TAD condition membership in {source}")
                if fraction is None or not math.isclose(fraction, support/n, rel_tol=1e-10, abs_tol=1e-12):
                    raise ValueError(f"Invalid TAD support fraction in {source}")
    ctx.notes.append(
        "TAD condition tracks reuse the workflow's recorded supported domains, including its boundary "
        "tolerance and replicate-support threshold. Their coordinates are representative observed intervals, "
        "not averaged coordinates. Saved condition strengths and mean insulation average the supporting "
        "called-domain observations, whose IDs and support/n_samples are retained in the source tables; "
        "they need not have complete coverage in every condition member. Unsupported tables combine insufficient support and overlapping "
        "segmentation conflicts; the upstream files do not record which reason applies to each row. "
        "Exact left/right boundary strengths, mean insulation and support remain in *.source.tsv. "
        "TAD blue=accepted, gray=unsupported or excluded. No significance is encoded.")


def _ada_color(value):
    if value is None:
        return "128,128,128"
    if value > 1:
        return "0,114,178"
    if value < 1:
        return "230,159,0"
    return "90,90,90"


def _export_ada(ctx):
    universe_source = "03_features/boundary_derived_domains/consensus.tsv"
    strength_source = "04_quantification/domain_scores/strength.tsv"
    if not ctx.exists(strength_source):
        if ctx.exists("05_aggregates/ada/aggregation.json"):
            raise ValueError("ADA aggregate metadata exists but domain strength table is missing")
        ctx.add_missing("ADA domain tracks", "The ADA domain-strength module is not available in this run copy.")
        return
    all_domains = _required(ctx, universe_source, ADA_FIELDS)
    _unique_ids(all_domains, universe_source)
    primary = int(ctx.config["INSULATION_PRIMARY_WINDOW_BP"])
    resolution = int(ctx.config["ADA_RESOLUTION"])
    domains = [row for row in all_domains if _integer(row["window_bp"], "window_bp") == primary]
    by_id = _unique_ids(domains, universe_source)
    quantitative_tracks = _nonoverlapping(domains)
    other = [dict(row, selection_status="outside_ADA_primary_window", reason="not_primary_insulation_window")
             for row in all_domains if _integer(row["window_bp"], "window_bp") != primary]
    ctx.copy_source(universe_source, "ada/source/all_boundary_derived_domains.tsv")
    ctx.write_table("ada/shared/attempted_domain_universe.tsv", domains,
                    fields=list(all_domains[0]) if all_domains else list(ADA_FIELDS))
    ctx.write_table("ada/shared/domains_outside_ADA_primary_window.tsv", other,
                    fields=(list(all_domains[0]) if all_domains else list(ADA_FIELDS)) + ["selection_status", "reason"])
    ctx.bed(rel="ada/shared/attempted_domain_universe.bed",
            rows=[(*_coordinate(row), row["feature_id"]) for row in domains],
            description=f"Exact ADA attempted domain universe shared by all observations; primary window {primary} bp",
            analysis="ada", scope="common", label="shared ADA domains", source=universe_source,
            resolution=resolution, window=primary)
    strengths = _required(ctx, strength_source, STRENGTH_FIELDS)
    by_sample = {sample["sample_id"]: {} for sample in ctx.samples}
    for row in strengths:
        sid, fid = row["sample_id"], row["feature_id"]
        if sid not in by_sample or fid not in by_id or fid in by_sample[sid]:
            raise ValueError(f"Unexpected or duplicate ADA sample/domain record: {sid}/{fid}")
        count = _integer(row["eligible_pixel_count"], "eligible_pixel_count")
        mean, total = _number(row["mean_intra_domain_oe"]), _number(row["oe_sum"])
        if count < 0 or total is None or total < 0:
            raise ValueError(f"Invalid ADA pixel count or O/E sum: {sid}/{fid}")
        if count == 0 and (mean is not None or total != 0):
            raise ValueError(f"ADA zero-coverage record contains a finite score or nonzero sum: {sid}/{fid}")
        if count > 0 and (mean is None or mean < 0 or not math.isclose(mean, total/count, rel_tol=1e-9, abs_tol=1e-12)):
            raise ValueError(f"ADA score and saved O/E sum disagree: {sid}/{fid}")
        by_sample[sid][fid] = row
    for sid, records in by_sample.items():
        if set(records) != set(by_id):
            raise ValueError(f"ADA domain universe is incomplete for {sid}; missing {len(set(by_id)-set(records))} records")
    ctx.copy_source(strength_source, "ada/source/domain_strengths.tsv")
    if ctx.exists("05_aggregates/ada/aggregation.json"):
        aggregation = ctx.json("05_aggregates/ada/aggregation.json")
        if int(aggregation.get("resolution", resolution)) != resolution:
            raise ValueError("ADA aggregate resolution differs from config")
        recorded_members = aggregation.get("condition_members", {})
        if recorded_members and {k: set(v) for k, v in recorded_members.items()} != {k: set(v) for k, v in ctx.conditions.items()}:
            raise ValueError("ADA aggregate condition membership differs from sample metadata")
        ctx.copy_source("05_aggregates/ada/aggregation.json", "ada/source/aggregation.json")
    row_fields = list(ADA_FIELDS) + ["sample_id", "mean_intra_domain_oe", "oe_sum", "eligible_pixel_count",
                                  "score_status", "score_missing_reason", "pileup_membership", "inference_status"]
    for sid, records in by_sample.items():
        rows, bed = [], []
        for domain in domains:
            fid = domain["feature_id"]
            measured = records[fid]
            score = _number(measured["mean_intra_domain_oe"])
            rows.append({**{key: domain[key] for key in ADA_FIELDS}, **measured,
                         "score_status": "finite_intra_domain_score" if score is not None else "unavailable_intra_domain_score",
                         "score_missing_reason": "" if score is not None else "no_eligible_upper_triangle_pixels",
                         "pileup_membership": "per_domain_resized_pixel_membership_not_recorded"})
            bed.append(_bed9(domain, _ada_color(score)))
        output = f"ada/individual/{_slug(sid)}/domain_strength"
        ctx.write_table(output + ".tsv", rows, fields=row_fields)
        ctx.bed(rel=output + ".bed", rows=bed,
                description=f"{sid}: exact ADA domains; O/E class colors; numerical values in matching TSV",
                analysis="ada", scope="individual", label=sid, source=strength_source,
                resolution=resolution, window=primary)
        if quantitative_tracks:
            ctx.bedgraph(rel=output + ".mean_intra_domain_oe.bedgraph",
                         rows=[(*_coordinate(row), float(row["mean_intra_domain_oe"]))
                               for row in rows if _number(row["mean_intra_domain_oe"]) is not None],
                         description=f"{sid}: saved mean intra-domain observed/expected; unavailable scores omitted",
                         analysis="ada", scope="individual", label=sid, source=strength_source,
                         resolution=resolution, window=primary)
            ctx.bedgraph(rel=output + ".eligible_pixel_count.bedgraph",
                         rows=[(*_coordinate(row), _integer(row["eligible_pixel_count"], "eligible_pixel_count"))
                               for row in rows],
                         description=f"{sid}: eligible upper-triangle pixels for each saved ADA strength; zero retained",
                         analysis="ada", scope="individual", label=sid, source=strength_source,
                         resolution=resolution, window=primary)
    summary_fields = list(ADA_FIELDS) + ["condition", "complete_case_mean_oe", "available_mean_oe", "available_sd_oe",
                                      "n_finite", "n_expected", "sample_ids", "missing_sample_ids", "mean_status",
                                      "score_missing_reason", "pileup_membership", "inference_status"]
    for condition, members in ctx.conditions.items():
        summaries, bed = [], []
        for domain in domains:
            fid = domain["feature_id"]
            values = {sid: _number(by_sample[sid][fid]["mean_intra_domain_oe"]) for sid in members}
            finite = [v for v in values.values() if v is not None]
            missing = [sid for sid, value in values.items() if value is None]
            available = statistics.mean(finite) if finite else None
            complete = available if not missing else None
            summaries.append({**{key: domain[key] for key in ADA_FIELDS}, "condition": condition,
                              "complete_case_mean_oe": complete if complete is not None else "NA",
                              "available_mean_oe": available if available is not None else "NA",
                              "available_sd_oe": statistics.stdev(finite) if len(finite) > 1 else "NA",
                              "n_finite": len(finite), "n_expected": len(members),
                              "sample_ids": ";".join(members), "missing_sample_ids": ";".join(missing),
                              "mean_status": "complete_case_mean" if not missing else "incomplete_case_mean_unavailable",
                              "score_missing_reason": "" if not missing else "one_or_more_replicates_have_no_eligible_upper_triangle_pixels",
                              "pileup_membership": "per_domain_resized_pixel_membership_not_recorded",
                              "inference_status": "descriptive"})
            bed.append(_bed9(domain, _ada_color(complete)))
        output = f"ada/conditions/{_slug(condition)}/domain_strength"
        ctx.write_table(output + ".tsv", summaries, fields=summary_fields)
        ctx.bed(rel=output + ".bed", rows=bed,
                description=f"{condition}: exact ADA domains; equal replicate mean only with complete score coverage; values in TSV",
                analysis="ada", scope="condition", label=condition, source=strength_source,
                resolution=resolution, window=primary)
        if quantitative_tracks:
            ctx.bedgraph(rel=output + ".complete_case_mean_oe.bedgraph",
                         rows=[(*_coordinate(row), float(row["complete_case_mean_oe"]))
                               for row in summaries if _number(row["complete_case_mean_oe"]) is not None],
                         description=f"{condition}: equal-replicate ADA strength mean where every observation has a finite score",
                         analysis="ada", scope="condition", label=condition, source=strength_source,
                         resolution=resolution, window=primary)
            ctx.bedgraph(rel=output + ".n_finite.bedgraph",
                         rows=[(*_coordinate(row), row["n_finite"]) for row in summaries],
                         description=f"{condition}: number of observations with finite intra-domain O/E; zero retained",
                         analysis="ada", scope="condition", label=condition, source=strength_source,
                         resolution=resolution, window=primary)
    if quantitative_tracks:
        ctx.notes.append(
            "The exact ADA primary-window domain set has no overlapping intervals, so numeric browser "
            "signals are also provided: individual mean intra-domain O/E and eligible-pixel counts, plus "
            "condition complete-case means and finite-observation counts. A domain's scalar value is "
            "displayed across that whole interval; it is not a per-bin measurement. Missing means are "
            "omitted, while genuine zero means and zero coverage counts are retained.")
    else:
        ctx.notes.append(
            "The exact ADA primary-window domain set contains overlapping intervals. Quantitative "
            "bedGraph tracks are therefore omitted without changing or merging the domains. The "
            "BED tracks retain every stable ID and exact interval, and the accompanying TSVs retain "
            "all numerical scores and coverage counts.")
    ctx.notes.append(
        "ADA attempted domains are exactly the primary-insulation-window rows of the saved "
        "boundary_derived_domains/consensus.tsv, evaluated for every original observation. They are adjacent "
        "all-sample consensus-boundary intervals and are distinct from the separate per-sample/condition "
        "TAD catalogue. Other window rows were not selected for this ADA. Strength is saved mean intra-domain "
        "observed/expected over eligible upper-triangle pixels. ADA BED color: blue O/E>1, orange O/E<1, "
        "dark gray O/E=1, gray unavailable. Numerical scores remain in matching TSVs; all exact intervals "
        "are retained, including unavailable scores. Condition colors use equal-observation means only when "
        "every member has a finite score. Sidecars also give available-member means, sample standard deviation "
        "and contribution counts; missing values are never zero-filled. BED always preserves exact "
        "domains; additional numeric bedGraph displays are created only for nonoverlapping domain sets. "
        "The resized ADA aggregate has pixel-specific "
        "valid_domain_count, not a scalar domain count; per-domain resized-pixel participation and detailed "
        "ADA exclusion reasons were not saved. A missing upper-triangle strength is not claimed to prove "
        "absence from every resized aggregate pixel.")
    ctx.add_missing("ADA per-domain resized-pixel membership and exclusion audit",
                    "Not recorded by the completed workflow. The exact attempted universe and saved score-coverage reasons are exported; no matrix queries or pileups were rerun.")


def export_domains(ctx):
    """Export existing TADs and the exact ADA candidate/score tables."""
    _export_tads(ctx)
    _export_ada(ctx)
