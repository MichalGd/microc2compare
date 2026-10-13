"""Insulation-defined domains and descriptive, replicate-supported changes.

Coordinates are zero-based half-open boundary-to-boundary intervals. Event
classes describe detectable geometry; they do not provide biological P/FDR.
"""
from __future__ import annotations

import math
import multiprocessing
import time
from bisect import bisect_left, bisect_right
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from .contracts import allocation_limits, atomic_json, digest, require, validate_metadata
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


class InsulationIndex:
    """View-local sorted arrays reused across thresholds and event queries."""

    def __init__(self, table, mask, c):
        require(len(table) == len(mask), "TAD score grid and registered mask differ")
        self.minimum = c["INSULATION_MIN_VALID_PIXEL_FRACTION"]
        self.views = {}
        self.windows = {}
        for key, frame in table.groupby(["chrom", "region"], sort=True):
            frame = frame.sort_values("start")
            valid = np.asarray(mask, dtype=bool)[frame.index]
            self.views[key] = (frame, frame.start.to_numpy(), valid,
                               np.r_[0, np.cumsum(~valid)])

    def window(self, key, window):
        cache_key = (key, window)
        if cache_key not in self.windows:
            frame = self.views[key][0]
            score = frame[f"log2_insulation_score_{window}"].to_numpy()
            good = np.isfinite(score) & (frame[f"valid_pixel_fraction_{window}"].to_numpy() >= self.minimum)
            self.windows[cache_key] = (score, good, np.r_[0, np.cumsum(good)])
        return self.windows[cache_key]

    def evaluable(self, key, start, end, window):
        if key not in self.views:
            return False
        _, starts, _, bad_mask = self.views[key]
        _, good, covered = self.window(key, window)
        a, b = np.searchsorted(starts, [start, end], side="left")
        ar, br = np.searchsorted(starts, [start, end], side="right")
        return bool(br > a and (ar-a)+(br-b) == 2 and
                    good[a:ar].all() and good[b:br].all() and
                    bad_mask[br] == bad_mask[a] and
                    (covered[br]-covered[a])/(br-a) >= self.minimum)


class IntervalIndex:
    """Existence of half-open overlaps, including nested rejected intervals."""

    def __init__(self, frames):
        self.views = {}
        frame = pd.concat(frames, ignore_index=True)
        for key, group in frame.groupby(["chrom", "region"], sort=False):
            group = group.sort_values("start")
            self.views[key] = (group.start.to_numpy(), np.maximum.accumulate(group.end.to_numpy()))

    def overlaps(self, key, start, end):
        if key not in self.views:
            return False
        starts, ends = self.views[key]
        last = int(np.searchsorted(starts, end, side="left"))
        return bool(last and ends[last-1] > start)


def call_domains(table, mask, c, window, sample_id, threshold=None, return_exclusions=False, *, index=None):
    """Call adjacent strong boundaries without crossing a excluded grid bin."""
    threshold = c["INSULATION_BOUNDARY_THRESHOLD"] if threshold is None else threshold
    score, strength = f"log2_insulation_score_{window}", f"boundary_strength_{window}"
    index = index or InsulationIndex(table, mask, c)
    rows, exclusions = [], []
    for (chrom, region), (frame, starts, valid, bad_mask) in index.views.items():
        scores, good, covered = index.window((chrom, region), window)
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
            a, b = np.searchsorted(starts, [start, end], side="left")
            ar, br = np.searchsorted(starts, [start, end], side="right")
            if b == a or bad_mask[b] != bad_mask[a]:
                exclude("excluded_grid_bins")
                continue
            if (ar-a)+(br-b)!=2 or not valid[a:ar].all() or not valid[b:br].all() or not good[a:ar].all() or not good[b:br].all():
                exclude("ineligible_boundary_bins")
                continue
            if (covered[b]-covered[a])/(b-a) < c["INSULATION_MIN_VALID_PIXEL_FRACTION"]:
                exclude("insufficient_score_coverage")
                continue
            rows.append({"feature_id": feature_id("tad", c["REFERENCE_ID"], (chrom, region, start, end, window), c["INSULATION_RESOLUTION"]),
                         "chrom": chrom, "start": start, "end": end, "region": region, "window_bp": window,
                         "length_bp": end-start, "left_boundary_strength": float(left[strength]),
                         "right_boundary_strength": float(right[strength]), "mean_log2_insulation": float(scores[a:b][good[a:b]].mean()),
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


def _sample_calls(task):
    sid, table, mask, c, windows, thresholds = task
    index = InsulationIndex(table, mask, c)
    return { (window, threshold): call_domains(table, mask, c, window, sid, threshold,
                                             return_exclusions=True, index=index)
             for window in windows for threshold in thresholds }


def worker_count(c, tables):
    """Admission estimate for resident parent tables and spawned worker copies.

    The estimate limits concurrency; it is not an OS-enforced memory limit.
    """
    limits = allocation_limits()
    sizes = [int(t.memory_usage(index=True, deep=True).sum()) for t in tables.values()]
    parent = sum(sizes)*4
    per_worker = max(max(sizes)*4, 64_000_000)
    memory = min(c.get("MEMORY_MB_TADS", 16000), c["MAX_MEMORY_MB"])*1_000_000
    if limits["available_memory_bytes"] is not None:
        memory = min(memory, limits["available_memory_bytes"])
    require(memory >= parent+per_worker, "TAD tables exceed the memory admission estimate; increase MEMORY_MB_TADS within the project allocation")
    return max(1, min(c.get("THREADS_TADS", 1), c["TOTAL_CPUS"],
                      max(1, int(limits["effective_cpus"])), len(tables),
                      (memory-parent)//per_worker))


def run(c, scope):
    from .science import save
    root = Path(c["run_dir"]); directory = root/"03_features/tads"
    samples = validate_metadata(c)["samples"]
    tables = {s["sample_id"]: pd.read_csv(root/f"03_features/insulation/{s['sample_id']}.tsv", sep="\t") for s in samples}
    with np.load(root/f"00_metadata/masks/{c['INSULATION_RESOLUTION']}.npz") as archive:
        mask = archive["mask"].copy()
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
    workers = worker_count(c, tables)
    started = time.perf_counter()
    print(f"TAD domain calling: {workers} worker(s), {len(tables)} samples, {len(windows)} windows, {len(thresholds)} thresholds", flush=True)
    tasks = [(sid, table, mask, c, windows, thresholds) for sid, table in tables.items()]
    if workers == 1:
        results = list(map(_sample_calls, tasks))
    else:
        # Spawn avoids inherited HDF5 handles and nested forked native libraries.
        with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn")) as pool:
            results = list(pool.map(_sample_calls, tasks))
    sample_calls = dict(zip(tables, results))
    domain_seconds = time.perf_counter()-started
    indices = {sid: InsulationIndex(table, mask, c) for sid, table in tables.items()}
    members = {cond: [s["sample_id"] for s in samples if s["condition"] == cond] for cond in conditions}

    def bed(frame, path, name_column):
        path.parent.mkdir(parents=True, exist_ok=True)
        frame[["chrom", "start", "end", name_column]].to_csv(path, sep="\t", header=False, index=False)
        if not len(frame): path.write_text("# No eligible domains or events\n")
        return path

    for window in windows:
        for threshold in thresholds:
            primary = threshold == primary_threshold
            calls={sid:result[(window,threshold)] for sid,result in sample_calls.items()}
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
                comparison_conditions = (comparison["reference_condition"], comparison["case_condition"])
                overlap_index = IntervalIndex([rejected[cond] for cond in comparison_conditions] +
                                              [calls[sid][1] for cond in comparison_conditions for sid in members[cond]])
                def evaluable(side, chrom, region, start, end):
                    ids = members[comparison[side+"_condition"]]
                    return sum(indices[sid].evaluable((chrom, region), start, end, window) for sid in ids)/len(ids)
                def ambiguous(chrom, region, start, end):
                    return overlap_index.overlaps((chrom, region), start, end)
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
    performance = directory/"performance.json"
    atomic_json(performance, {"workers": workers, "samples": len(tables),
                             "windows": len(windows), "thresholds": len(thresholds),
                             "domain_calling_seconds": domain_seconds,
                             "total_seconds": time.perf_counter()-started,
                             "memory_admission": "conservative estimate; not an OS limit"})
    outputs.append(performance)
    return outputs
