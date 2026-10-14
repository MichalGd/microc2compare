"""Insulation-defined domains and descriptive, replicate-supported changes.

Coordinates are zero-based half-open boundary-to-boundary intervals. Event
classes describe detectable geometry; they do not provide biological P/FDR.
"""
from __future__ import annotations

import math
from bisect import bisect_left, bisect_right
from pathlib import Path

import numpy as np
import pandas as pd

from .contracts import atomic_json, digest, require, validate_metadata
from .features import feature_id

DOMAIN_COLUMNS = ["feature_id", "chrom", "start", "end", "region", "window_bp", "length_bp",
                  "left_boundary_strength", "right_boundary_strength", "mean_log2_insulation",
                  "support", "n_samples", "support_fraction", "sample_ids", "inference_status"]
EVENT_COLUMNS = ["event_id", "event", "comparison_id", "reference_condition", "case_condition",
                 "chrom", "start", "end", "region", "window_bp", "reference_ids", "case_ids",
                 "reference_start", "reference_end", "case_start", "case_end", "reference_length_bp",
                 "case_length_bp", "length_change_bp", "length_ratio", "start_shift_bp", "end_shift_bp",
                 "reference_support_fraction", "case_support_fraction", "reference_evaluable_fraction",
                 "case_evaluable_fraction", "boundary_strength_change", "mean_log2_insulation_change",
                 "reciprocal_overlap", "reason", "inference_status", "PValue", "FDR"]


def required_support(n, fraction):
    return max(1, math.ceil(n * fraction))


def call_domains(table, mask, c, window, sample_id, threshold=None, return_exclusions=False):
    """Call adjacent strong boundaries without crossing a excluded grid bin."""
    threshold = c["INSULATION_BOUNDARY_THRESHOLD"] if threshold is None else threshold
    score, strength = f"log2_insulation_score_{window}", f"boundary_strength_{window}"
    require(len(table) == len(mask), "TAD score grid and registered mask differ")
    rows, exclusions = [], []
    for (chrom, region), frame in table.groupby(["chrom", "region"], sort=True):
        frame = frame.sort_values("start")
        if threshold in ("Li", "Otsu"):
            require(threshold == c["INSULATION_BOUNDARY_THRESHOLD"], "Automatic threshold sensitivity needs separate insulation computation")
            called = frame[f"is_boundary_{window}"].fillna(False).astype(bool)
        else: called = frame[strength] >= float(threshold)
        boundaries = frame.loc[np.isfinite(frame[strength]) & called]
        for left, right in zip(boundaries.itertuples(index=False, name=None), boundaries.iloc[1:].itertuples(index=False, name=None)):
            left, right = dict(zip(boundaries.columns, left)), dict(zip(boundaries.columns, right))
            start, end = int(left["start"]), int(right["start"])
            def exclude(reason):
                exclusions.append(dict(chrom=chrom,start=start,end=end,region=region,window_bp=window,reason=reason))
            if not c["TAD_MIN_LENGTH_BP"] <= end - start <= c["TAD_MAX_LENGTH_BP"]:
                exclude("outside_domain_length_limits")
                continue
            inside = frame.loc[(frame.start >= start) & (frame.start < end)]
            if not len(inside) or not np.asarray(mask)[inside.index].all():
                exclude("excluded_grid_bins")
                continue
            outer=frame.loc[frame.start.isin([start,end])]
            if len(outer)!=2 or not np.asarray(mask)[outer.index].all() or not np.isfinite(outer[score]).all() or not (outer[f"valid_pixel_fraction_{window}"]>=c["INSULATION_MIN_VALID_PIXEL_FRACTION"]).all():
                exclude("ineligible_boundary_bins")
                continue
            good = np.isfinite(inside[score]) & (inside[f"valid_pixel_fraction_{window}"] >= c["INSULATION_MIN_VALID_PIXEL_FRACTION"])
            if good.mean() < c["INSULATION_MIN_VALID_PIXEL_FRACTION"]:
                exclude("insufficient_score_coverage")
                continue
            rows.append({"feature_id": feature_id("tad", c["REFERENCE_ID"], (chrom, region, start, end, window), c["INSULATION_RESOLUTION"]),
                         "chrom": chrom, "start": start, "end": end, "region": region, "window_bp": window,
                         "length_bp": end-start, "left_boundary_strength": float(left[strength]),
                         "right_boundary_strength": float(right[strength]), "mean_log2_insulation": float(inside.loc[good, score].mean()),
                         "support": 1, "n_samples": 1, "support_fraction": 1.0,
                         "sample_ids": sample_id, "inference_status": "descriptive"})
    domains=pd.DataFrame(rows,columns=DOMAIN_COLUMNS)
    rejected=pd.DataFrame(exclusions,columns=["chrom","start","end","region","window_bp","reason"])
    return (domains,rejected) if return_exclusions else domains


def condition_domains(frames, c, window):
    """Complete-span coordinate clusters; only distinct biological samples vote."""
    points = sorted((row.chrom, row.region, int(row.start), int(row.end), sid, row)
                    for sid, frame in frames.items() for row in frame.itertuples())
    groups, active, current_view = [], [], None
    tolerance = c["TAD_BOUNDARY_TOLERANCE_BP"]
    for point in points:
        if point[:2] != current_view: active = []; current_view = point[:2]
        active = [g for g in active if point[2]-g[0][2] <= tolerance]
        matches = [g for g in active if g[0][:2] == point[:2]
                   and point[4] not in {p[4] for p in g}
                   and max([p[2] for p in g] + [point[2]]) - min([p[2] for p in g] + [point[2]]) <= tolerance
                   and max([p[3] for p in g] + [point[3]]) - min([p[3] for p in g] + [point[3]]) <= tolerance]
        if matches:
            min(matches, key=lambda g: (sum(abs(point[2]-p[2])+abs(point[3]-p[3]) for p in g), g[0][:5])).append(point)
        else:
            group=[point]; groups.append(group); active.append(group)
    accepted, rejected = [], []
    for group in groups:
        representative = min(group, key=lambda p: (sum(abs(p[2]-q[2])+abs(p[3]-q[3]) for q in group), p[:5]))
        chrom, region, start, end = representative[:4]
        strengths = [p[5] for p in group]
        row = {"feature_id": feature_id("tad", c["REFERENCE_ID"], (chrom, region, start, end, window), c["INSULATION_RESOLUTION"]),
               "chrom": chrom, "start": start, "end": end, "region": region, "window_bp": window,
               "length_bp": end-start, "left_boundary_strength": float(np.mean([p.left_boundary_strength for p in strengths])),
               "right_boundary_strength": float(np.mean([p.right_boundary_strength for p in strengths])),
               "mean_log2_insulation": float(np.mean([p.mean_log2_insulation for p in strengths])),
               "support": len(group), "n_samples": len(frames), "support_fraction": len(group)/len(frames),
               "sample_ids": ";".join(sorted(p[4] for p in group)), "inference_status": "descriptive"}
        (accepted if len(group) >= required_support(len(frames), c["TAD_MIN_REPLICATE_FRACTION"]) else rejected).append(row)
    # Overlapping supported condition intervals imply inconsistent segmentation.
    # Retain them in the disagreement table rather than hiding the ambiguity.
    conflicted = set()
    active=[]; current_view=None
    for i, a in sorted(enumerate(accepted), key=lambda item:(item[1]["chrom"],item[1]["region"],item[1]["start"])):
        key=(a["chrom"],a["region"])
        if key!=current_view:active=[];current_view=key
        active=[(j,b) for j,b in active if b["end"]>a["start"]]
        for j,b in active:conflicted.update((i,j))
        active.append((i,a))
    rejected += [r for i, r in enumerate(accepted) if i in conflicted]
    accepted = [r for i, r in enumerate(accepted) if i not in conflicted]
    return pd.DataFrame(accepted, columns=DOMAIN_COLUMNS), pd.DataFrame(rejected, columns=DOMAIN_COLUMNS)


def overlap(a, b):
    return max(0, min(a["end"], b["end"])-max(a["start"], b["start"]))


def match_domains(reference, case, c, comparison, window, evaluable=None, ambiguous=None):
    """Classify connected correspondences once, including coherent split/merge.

    Small-over-large overlap admits split/merge edges. One-to-one matches then
    require reciprocal overlap. Ambiguous groups never become disappearance.
    """
    left, right = reference.to_dict("records"), case.to_dict("records")
    nodes = [("r", i) for i in range(len(left))] + [("c", i) for i in range(len(right))]
    links = {n: set() for n in nodes}
    by_view={}
    for j,b in enumerate(right):by_view.setdefault((b["chrom"],b["region"]),[]).append((j,b))
    indices={}
    for key,items in by_view.items():
        items.sort(key=lambda pair:pair[1]["start"])
        starts=[b["start"] for j,b in items];ends=[b["end"] for j,b in items]
        indices[key]=(items,starts,ends,ends==sorted(ends))
    for i, a in enumerate(left):
        if (a["chrom"],a["region"]) not in indices:continue
        items,starts,ends,ordered_ends=indices[(a["chrom"],a["region"])]
        first=bisect_right(ends,a["start"]) if ordered_ends else 0
        last=bisect_left(starts,a["end"])
        for j, b in items[first:last]:
            shared = overlap(a, b)
            if shared/min(a["end"]-a["start"], b["end"]-b["start"]) >= c["TAD_MIN_RECIPROCAL_OVERLAP"]:
                links[("r", i)].add(("c", j)); links[("c", j)].add(("r", i))
    groups, visited = [], set()
    for node in nodes:
        if node in visited: continue
        pending, group = [node], []
        while pending:
            item = pending.pop()
            if item in visited: continue
            visited.add(item); group.append(item); pending.extend(sorted(links[item]))
        groups.append(group)
    rows = []
    tolerance = c["TAD_BOUNDARY_TOLERANCE_BP"]
    def contiguous(items):
        ordered = sorted(items, key=lambda x: x["start"])
        return all(abs(a["end"]-b["start"]) <= tolerance for a, b in zip(ordered, ordered[1:]))
    for group in groups:
        r = [left[i] for side, i in group if side == "r"]
        t = [right[i] for side, i in group if side == "c"]
        items = r + t; chrom, region = items[0]["chrom"], items[0]["region"]
        start, end = min(x["start"] for x in items), max(x["end"] for x in items)
        rf = evaluable("reference", chrom, region, start, end) if evaluable else 1.0
        cf = evaluable("case", chrom, region, start, end) if evaluable else 1.0
        event, reason, reciprocal = "complex", "ambiguous_correspondence", np.nan
        if min(rf, cf) < c["TAD_MIN_REPLICATE_FRACTION"]:
            event, reason = "uncertain", "insufficient_common_coverage_or_view_edge"
        elif ambiguous and ambiguous(chrom, region, start, end):
            event, reason = "uncertain", "replicate_domain_disagreement"
        elif not r:
            event, reason = "appearing", "newly_detected_not_explained_by_split"
        elif not t:
            event, reason = "disappearing", "no_longer_detected_not_explained_by_merge"
        elif len(r) == 1 and len(t) == 1:
            a, b = r[0], t[0]; la, lb = a["end"]-a["start"], b["end"]-b["start"]
            reciprocal = min(overlap(a, b)/la, overlap(a, b)/lb)
            if reciprocal >= c["TAD_MIN_RECIPROCAL_OVERLAP"]:
                delta = lb-la
                if abs(delta) > c["TAD_SIZE_CHANGE_BP"] and abs(delta)/la > c["TAD_SIZE_CHANGE_FRACTION"]:
                    event = "growing" if delta > 0 else "shrinking"
                elif max(abs(b["start"]-a["start"]), abs(b["end"]-a["end"])) > tolerance:
                    event = "shifted"
                else: event = "stable"
                reason = "one_to_one_match"
        elif len(r) > 1 and len(t) == 1 or len(t) > 1 and len(r) == 1:
            many, one = (r, t[0]) if len(r) > 1 else (t, r[0])
            outer = abs(min(x["start"] for x in many)-one["start"]) <= tolerance and abs(max(x["end"] for x in many)-one["end"]) <= tolerance
            covered = sum(overlap(x, one) for x in many)/(one["end"]-one["start"])
            if contiguous(many) and outer and covered >= c["TAD_MIN_RECIPROCAL_OVERLAP"]:
                event = "merge" if len(r) > 1 else "split"
                reason = "coherent_adjacent_domain_correspondence"
        ra, rb = (min(x["start"] for x in r), max(x["end"] for x in r)) if r else (np.nan, np.nan)
        ca, cb = (min(x["start"] for x in t), max(x["end"] for x in t)) if t else (np.nan, np.nan)
        rs, cs = rb-ra, cb-ca
        strength = lambda xs: float(np.mean([(x["left_boundary_strength"]+x["right_boundary_strength"])/2 for x in xs])) if xs else np.nan
        strength_change = strength(t)-strength(r)
        if event == "stable" and np.isfinite(strength_change) and abs(strength_change) > c["TAD_STRENGTH_CHANGE"]:
            event = "strength_change"
        info = {"chrom": chrom, "start": start, "end": end, "region": region, "window_bp": window,
                "reference_ids": ";".join(sorted(x["feature_id"] for x in r)), "case_ids": ";".join(sorted(x["feature_id"] for x in t))}
        row = dict(info, event_id="event_"+digest([comparison["comparison_id"], info])[:16], event=event,
                   comparison_id=comparison["comparison_id"], reference_condition=comparison["reference_condition"], case_condition=comparison["case_condition"],
                   reference_start=ra, reference_end=rb, case_start=ca, case_end=cb, reference_length_bp=rs, case_length_bp=cs,
                   length_change_bp=cs-rs, length_ratio=cs/rs if np.isfinite(rs) else np.nan,
                   start_shift_bp=ca-ra, end_shift_bp=cb-rb,
                   reference_support_fraction=min(x["support_fraction"] for x in r) if r else np.nan,
                   case_support_fraction=min(x["support_fraction"] for x in t) if t else np.nan,
                   reference_evaluable_fraction=rf, case_evaluable_fraction=cf, boundary_strength_change=strength_change,
                   mean_log2_insulation_change=(float(np.mean([x["mean_log2_insulation"] for x in t]))-float(np.mean([x["mean_log2_insulation"] for x in r]))) if r and t else np.nan,
                   reciprocal_overlap=reciprocal, reason=reason, inference_status="descriptive", PValue=np.nan, FDR=np.nan)
        rows.append(row)
    return pd.DataFrame(rows, columns=EVENT_COLUMNS)


def run(c, scope):
    from .science import save
    root = Path(c["run_dir"]); directory = root/"03_features/tads"
    samples = validate_metadata(c)["samples"]
    tables = {s["sample_id"]: pd.read_csv(root/f"03_features/insulation/{s['sample_id']}.tsv", sep="\t") for s in samples}
    archive = np.load(root/f"00_metadata/masks/{c['INSULATION_RESOLUTION']}.npz")
    bins = pd.read_csv(root/f"00_metadata/masks/{c['INSULATION_RESOLUTION']}.bins.tsv", sep="\t")
    for sid, table in tables.items():
        require(table[["chrom", "start", "end"]].equals(bins[["chrom", "start", "end"]]), f"{sid}: TAD score grid differs from registered bins")
    comparisons = validate_metadata(c)["comparisons"]
    conditions = sorted({s["condition"] for s in samples})
    condition_ids = {cond: "condition_"+digest(cond)[:12] for cond in conditions}
    outputs = []; sensitivity = []
    windows = [int(w) for w in c["INSULATION_WINDOWS_BP"].split(",")]
    primary_threshold = c["INSULATION_BOUNDARY_THRESHOLD"]
    if primary_threshold not in ("Li", "Otsu"): primary_threshold = float(primary_threshold)
    thresholds = [primary_threshold] + sorted(set(float(t) for t in c["TAD_SENSITIVITY_THRESHOLDS"].split(",") if t) - ({primary_threshold} if isinstance(primary_threshold, float) else set()))

    def bed(frame, path, name_column):
        path.parent.mkdir(parents=True, exist_ok=True)
        frame[["chrom", "start", "end", name_column]].to_csv(path, sep="\t", header=False, index=False)
        if not len(frame): path.write_text("# No eligible domains or events\n")
        return path

    for window in windows:
        for threshold in thresholds:
            primary = threshold == primary_threshold
            calls={sid:call_domains(table,archive["mask"],c,window,sid,threshold,return_exclusions=True) for sid,table in tables.items()}
            per_sample={sid:result[0] for sid,result in calls.items()}
            supported, rejected = {}, {}
            for cond in conditions:
                ids = [s["sample_id"] for s in samples if s["condition"] == cond]
                supported[cond], rejected[cond] = condition_domains({sid: per_sample[sid] for sid in ids}, c, window)
                if primary:
                    base = directory/"conditions"/f"{condition_ids[cond]}.{window}"
                    outputs.extend([save(supported[cond], Path(str(base)+".tsv")), save(rejected[cond], Path(str(base)+".unsupported.tsv")), bed(supported[cond], Path(str(base)+".bed"), "feature_id")])
            if primary:
                for sid, frame in per_sample.items():
                    base = directory/"samples"/f"{sid}.{window}"
                    outputs.extend([save(frame, Path(str(base)+".tsv")), bed(frame, Path(str(base)+".bed"), "feature_id")])
                    outputs.append(save(calls[sid][1],Path(str(base)+".exclusions.tsv")))
            for comparison in comparisons:
                def evaluable(side, chrom, region, start, end):
                    cond = comparison[side+"_condition"]; ids = [s["sample_id"] for s in samples if s["condition"] == cond]
                    good_samples = 0
                    for sid in ids:
                        table = tables[sid]
                        area = table.loc[(table.chrom == chrom) & (table.region == region) & (table.start >= start) & (table.start <= end)]
                        score = f"log2_insulation_score_{window}"
                        ok = np.isfinite(area[score]) & (area[f"valid_pixel_fraction_{window}"] >= c["INSULATION_MIN_VALID_PIXEL_FRACTION"])
                        # Both outer boundary bins must be present and usable.
                        outer = area.loc[area.start.isin([start, end])]
                        if len(area) and len(outer) == 2 and np.isfinite(outer[score]).all() and (outer[f"valid_pixel_fraction_{window}"] >= c["INSULATION_MIN_VALID_PIXEL_FRACTION"]).all() and archive["mask"][area.index].all() and ok.mean() >= c["INSULATION_MIN_VALID_PIXEL_FRACTION"]:
                            good_samples += 1
                    return good_samples/len(ids)
                def ambiguous(chrom, region, start, end):
                    for cond in (comparison["reference_condition"], comparison["case_condition"]):
                        frame = rejected[cond]
                        if ((frame.chrom == chrom) & (frame.region == region) & (frame.start < end) & (frame.end > start)).any(): return True
                        for sample in samples:
                            if sample["condition"]!=cond:continue
                            excluded=calls[sample["sample_id"]][1]
                            if ((excluded.chrom==chrom)&(excluded.region==region)&(excluded.start<end)&(excluded.end>start)).any():return True
                    return False
                events = match_domains(supported[comparison["reference_condition"]], supported[comparison["case_condition"]], c, comparison, window, evaluable, ambiguous)
                for event, count in events.event.value_counts().items():
                    sensitivity.append({"comparison_id": comparison["comparison_id"], "window_bp": window, "boundary_threshold": threshold, "primary": primary, "event": event, "count": int(count), "inference_status": "descriptive"})
                if not len(events): sensitivity.append({"comparison_id": comparison["comparison_id"], "window_bp": window, "boundary_threshold": threshold, "primary": primary, "event": "no_events", "count": 0, "inference_status": "descriptive"})
                if primary:
                    base = root/"06_comparisons"/comparison["comparison_id"]/f"tad_events.{window}"
                    outputs.extend([save(events, Path(str(base)+".tsv")), bed(events, Path(str(base)+".bed"), "event_id")])
                    counts = events.event.value_counts().rename_axis("event").reset_index(name="count")
                    outputs.append(save(counts, Path(str(base)+".summary.tsv")))
                else:
                    outputs.append(save(events, directory/"sensitivity"/f"{comparison['comparison_id']}.{window}.{threshold:g}.events.tsv"))
    outputs.append(save(pd.DataFrame(sensitivity, columns=["comparison_id", "window_bp", "boundary_threshold", "primary", "event", "count", "inference_status"]), directory/"sensitivity.tsv"))
    parameters = directory/"parameters.json"
    atomic_json(parameters, {"method": "insulation_defined_adjacent_boundaries", "conditions": condition_ids,
                            "parameters": {k: v for k, v in c.items() if k.startswith(("TAD_", "INSULATION_"))},
                            "inference_status": "descriptive", "event_fdr": "not_available",
                            "coordinate_convention": "zero_based_half_open", "thresholds": thresholds})
    outputs.append(parameters)
    return outputs


def legacy_quantify(clr, features, bins, mask, labels, c):
    lookup = {(str(b.chrom), int(b.start)): int(i) for i, b in enumerate(bins.itertuples())}
    targets, reasons = {}, {}
    for index, row in features.iterrows():
        i = lookup.get((str(row.chrom1), int(row.start1)))
        j = lookup.get((str(row.chrom2), int(row.start2)))
        reason = None
        if i is None or j is None: reason = "outside_grid"
        elif i == j: reason = "collapsed_to_diagonal"
        elif row.chrom1 != row.chrom2: reason = "trans_unsupported"
        elif not mask[i] or not mask[j] or labels[i] < 0 or labels[i] != labels[j]: reason = "masked_or_different_view"
        elif not c["LOOP_MIN_DISTANCE_BP"] <= int(row.start2-row.start1) <= c["LOOP_MAX_DISTANCE_BP"]: reason = "outside_distance_range"
        if reason: reasons[index] = reason
        else: targets[(min(i, j), max(i, j))] = index
    counts = np.zeros(len(features), dtype=np.int64)
    exposure = 0
    positions = bins.start.to_numpy()
    for chunk in pixel_chunks(clr, c["PIXEL_CHUNKSIZE"]):
        i, j, y = [chunk[k].to_numpy() for k in ("bin1_id", "bin2_id", "count")]
        distances = positions[j]-positions[i]
        eligible = mask[i] & mask[j] & (labels[i] >= 0) & (labels[i] == labels[j]) & (distances >= c["LOOP_MIN_DISTANCE_BP"]) & (distances <= c["LOOP_MAX_DISTANCE_BP"])
        exposure += int(y[eligible].sum())
        for a, b, value in zip(i, j, y):
            index = targets.get((int(a), int(b)))
            if index is not None: counts[index] += int(value)
    require(exposure > 0, "Eligible cis-contact exposure is zero; cannot quantify relative contact allocation")
    values = pd.Series(counts, dtype="Int64")
    for index in reasons: values.loc[index] = pd.NA
    return values, exposure, reasons


def legacy_oe_window(clr, region1, region2, values, labels, views, expected, c):
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

from .matrix import pixel_chunks, bounded_raw


def legacy_consensus_loops(frames, reference, resolution, tolerance=1):
    records = []
    for sample, frame in frames.items():
        if len(frame):
            for row in canonical_loops(frame, reference, resolution).to_dict("records"):
                row["source"] = sample
                records.append(row)
    records.sort(key=lambda r: (r["chrom1"], r["chrom2"], r["start1"], r["start2"], r["source"]))
    groups = []
    for record in records:
        chosen = None
        for group in groups:
            if all(record["chrom1"] == x["chrom1"] and record["chrom2"] == x["chrom2"] and
                   abs(record["start1"]-x["start1"]) <= resolution*tolerance and
                   abs(record["start2"]-x["start2"]) <= resolution*tolerance for x in group):
                chosen = group; break
        if chosen is None: groups.append([record])
        else: chosen.append(record)
    rows = []
    for group in groups:
        rep = min(group, key=lambda a: (sum(abs(a["start1"]-b["start1"])+abs(a["start2"]-b["start2"]) for b in group), a["start1"], a["start2"], a["source"]))
        row = {key: rep[key] for key in LOOP_COLUMNS}
        biological = sorted({r["source"] for r in group if not r["source"].startswith("pool__")})
        row.update({"member_sources": ";".join(sorted({r["source"] for r in group})), "biological_support": len(biological), "inference_status": "descriptive_same_data"})
        rows.append(row)
    return pd.DataFrame(rows, columns=LOOP_COLUMNS+["member_sources", "biological_support", "inference_status"])
from .features import canonical_loops, LOOP_COLUMNS
