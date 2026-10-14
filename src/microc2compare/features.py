"""Deterministic candidate coordinates and exact raw-contact quantification."""
from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd

from .contracts import require
from .matrix import pixel_chunks

LOOP_COLUMNS = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "feature_id"]


def feature_id(kind, reference, coordinates, resolution):
    value = "|".join(map(str, [kind, reference, *coordinates, resolution]))
    return kind + "_" + hashlib.sha256(value.encode()).hexdigest()[:16]


def canonical_loops(frame, reference, resolution):
    required = LOOP_COLUMNS[:6]
    require(set(required) <= set(frame), "BEDPE requires chrom1,start1,end1,chrom2,start2,end2 headers")
    rows = []
    for row in frame.to_dict("records"):
        a = (str(row["chrom1"]), int(row["start1"]), int(row["end1"]))
        b = (str(row["chrom2"]), int(row["start2"]), int(row["end2"]))
        require(a[1] >= 0 and b[1] >= 0 and a[2] > a[1] and b[2] > b[1], "Invalid BEDPE coordinates")
        require(a[1] % resolution == b[1] % resolution == 0 and a[2]-a[1] == b[2]-b[1] == resolution,
                "v0.1 loop anchors must match the declared single-bin grid; no silent conversion")
        if b < a: a, b = b, a
        require(not (a[0] == b[0] and a[2] > b[1]), "Overlapping/diagonal loop anchors are unsupported")
        values = a + b
        row.update(dict(zip(required, values)))
        row["feature_id"] = feature_id("loop", reference, values, resolution)
        rows.append(row)
    result = pd.DataFrame(rows, columns=list(dict.fromkeys(list(frame.columns) + ["feature_id"])))
    return result.drop_duplicates(required).sort_values(required).reset_index(drop=True)


def consensus_loops(frames, reference, resolution, tolerance=1):
    records = []
    for sample, frame in frames.items():
        if len(frame):
            for row in canonical_loops(frame, reference, resolution).to_dict("records"):
                row["source"] = sample
                records.append(row)
    records.sort(key=lambda r: (r["chrom1"], r["chrom2"], r["start1"], r["start2"], r["source"]))
    groups, active = [], []
    for record in records:
        # Earlier clusters outside the first-anchor tolerance cannot match any
        # later sorted record. Keep original first-match complete-linkage order.
        active = [group for group in active
                  if group[0]["chrom1"] == record["chrom1"] and group[0]["chrom2"] == record["chrom2"]
                  and record["start1"]-group[0]["start1"] <= resolution*tolerance]
        chosen = None
        for group in active:
            if all(record["chrom1"] == x["chrom1"] and record["chrom2"] == x["chrom2"] and
                   abs(record["start1"]-x["start1"]) <= resolution*tolerance and
                   abs(record["start2"]-x["start2"]) <= resolution*tolerance for x in group):
                chosen = group; break
        if chosen is None:
            group = [record]; groups.append(group); active.append(group)
        else: chosen.append(record)
    rows = []
    for group in groups:
        rep = min(group, key=lambda a: (sum(abs(a["start1"]-b["start1"])+abs(a["start2"]-b["start2"]) for b in group), a["start1"], a["start2"], a["source"]))
        row = {key: rep[key] for key in LOOP_COLUMNS}
        biological = sorted({r["source"] for r in group if not r["source"].startswith("pool__")})
        row.update({"member_sources": ";".join(sorted({r["source"] for r in group})), "biological_support": len(biological), "inference_status": "descriptive_same_data"})
        rows.append(row)
    return pd.DataFrame(rows, columns=LOOP_COLUMNS+["member_sources", "biological_support", "inference_status"])


def quantify(clr, features, bins, mask, labels, c):
    return quantify_many(clr, {"features": features}, bins, mask, labels, c)["features"]


def quantify_many(clr, universes, bins, mask, labels, c):
    """Measure candidate sets in one bounded, vectorized raw-pixel pass.

    Exposure, zero counts, missing counts and exclusion rules are unchanged.
    Coordinates use int64 keys; no dense chromosome/genome matrix is allocated.
    """
    lookup = {(str(b.chrom), int(b.start)): int(i) for i, b in enumerate(bins.itertuples())}
    n = len(bins)
    require(n*n <= np.iinfo(np.int64).max, "Grid exceeds int64 pixel-key capacity")
    sets = {}
    for name, features in universes.items():
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
            else: targets[min(i, j)*n+max(i, j)] = index
        sets[name] = (targets, reasons)
    keys = np.array(sorted({key for targets, _ in sets.values() for key in targets}), dtype=np.int64)
    measured = np.zeros(len(keys), dtype=np.int64)
    exposure = 0
    positions = bins.start.to_numpy()
    for chunk in pixel_chunks(clr, c["PIXEL_CHUNKSIZE"]):
        i, j, y = [chunk[k].to_numpy() for k in ("bin1_id", "bin2_id", "count")]
        distances = positions[j]-positions[i]
        eligible = mask[i] & mask[j] & (labels[i] >= 0) & (labels[i] == labels[j]) & (distances >= c["LOOP_MIN_DISTANCE_BP"]) & (distances <= c["LOOP_MAX_DISTANCE_BP"])
        exposure += int(y[eligible].sum())
        if len(keys):
            pixel_keys = i.astype(np.int64)*n+j.astype(np.int64)
            # Restrict the search to the candidate key range, then exact matches.
            in_range = (pixel_keys >= keys[0]) & (pixel_keys <= keys[-1])
            selected, selected_y = pixel_keys[in_range], y[in_range]
            at = np.searchsorted(keys, selected)
            matched = keys[at] == selected
            np.add.at(measured, at[matched], selected_y[matched])
    require(exposure > 0, "Eligible cis-contact exposure is zero; cannot quantify relative contact allocation")
    result = {}
    for name, (targets, reasons) in sets.items():
        counts = np.zeros(len(universes[name]), dtype=np.int64)
        if targets:
            counts[list(targets.values())] = measured[np.searchsorted(keys, list(targets))]
        values = pd.Series(counts, dtype="Int64")
        for index in reasons: values.loc[index] = pd.NA
        result[name] = (values, exposure, reasons)
    return result


def boundary_consensus(frames, reference, resolution, window):
    points = []
    for sample, frame in frames.items():
        for row in frame.itertuples():
            points.append((row.chrom, row.region, int(row.start), sample))
    groups = []
    for point in sorted(points):
        if groups and point[:2] == groups[-1][0][:2] and point[2]-groups[-1][0][2] <= resolution:
            groups[-1].append(point)
        else: groups.append([point])
    rows = []
    for group in groups:
        chrom, view = group[0][:2]
        start = int(np.floor(np.median([p[2] for p in group])/resolution+0.5))*resolution
        rows.append({"chrom": chrom, "start": start, "end": start+resolution, "region": view,
                     "window_bp": window, "support": len({p[3] for p in group if not p[3].startswith("pool__")}),
                     "feature_id": feature_id("boundary", reference, (chrom, start, view, window), resolution)})
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "region", "window_bp", "support", "feature_id"])
