"""Optional, pinned CPU stripe discovery on contiguous common eligible views.

StripePy's native HDF5 bounds are inclusive bin indices. Canonical exports use
zero-based half-open genomic intervals, and retain those native bounds explicitly.
"""
from __future__ import annotations

import importlib.metadata
from contextlib import contextmanager
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import cooler
import numpy as np
import pandas as pd

from .contracts import atomic_json, digest, file_hash, read_json, require, validate_metadata
from .matrix import pixel_chunks

SUPPORTED_VERSION = "1.3.0"
BEDPE_COLUMNS = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
STRIPE_COLUMNS = BEDPE_COLUMNS + [
    "stripe_id", "sample_id", "condition", "reference_id", "resolution",
    "anchor_chrom", "anchor_start", "anchor_end", "body_chrom", "body_start", "body_end",
    "seed_bp", "direction", "matrix_location", "width_bp", "length_bp", "view", "segment_id",
    "native_left_bin", "native_right_bin", "native_top_bin", "native_bottom_bin",
    "relative_change", "top_persistence", "inner_mean", "outer_mean", "coefficient_of_variation",
    "inference_status",
]
SEGMENT_COLUMNS = ["segment_id", "chrom", "start", "end", "view", "native_first_bin", "native_last_bin", "bin_count", "status"]
EXCLUSION_COLUMNS = ["segment_id", "matrix_location", "seed_bin", "reason"]


def stripe_resolutions(c):
    parts = str(c.get("STRIPE_RESOLUTIONS", "4000")).split(",")
    require(parts and all(p.isdigit() and int(p) > 0 for p in parts), "STRIPE_RESOLUTIONS requires positive integer grids")
    require(len(set(parts)) == len(parts), "Duplicate stripe resolutions")
    return sorted(int(p) for p in parts)


def probe():
    try:
        version = importlib.metadata.version("stripepy-hic")
    except importlib.metadata.PackageNotFoundError as error:
        raise ValueError("Stripe detection requires the optional pinned stripe environment; install envs/stripes.yml") from error
    require(version == SUPPORTED_VERSION, f"Stripe adapter requires stripepy-hic=={SUPPORTED_VERSION}; found {version}")
    from stripepy.data_structures import ResultFile
    require(callable(ResultFile), "StripePy ResultFile API unavailable")
    return {"software": "stripepy-hic", "version": version, "normalization": "raw_counts", "cpu_only": True}


def segments_for_context(bins, mask, labels, views, resolution):
    """Never turn unavailable bins into zero signal: split at every masked gap."""
    rows = []
    for vid, view in enumerate(views.itertuples(index=False)):
        ids = np.flatnonzero((labels == vid) & mask)
        cuts = np.flatnonzero(np.diff(ids) != 1) + 1
        for group in np.split(ids, cuts):
            if not len(group):
                continue
            first, last = int(group[0]), int(group[-1])
            start, end = int(bins.iloc[first].start), int(bins.iloc[last].end)
            require((bins.iloc[group].chrom.astype(str) == str(view.chrom)).all(), "Stripe segment crosses chromosomes")
            rows.append(dict(segment_id=f"segment{len(rows):06d}", chrom=str(view.chrom), start=start, end=end,
                             view=str(view.name), native_first_bin=first, native_last_bin=last,
                             bin_count=len(group), status="eligible" if len(group) >= 3 else "too_short"))
    return pd.DataFrame(rows, columns=SEGMENT_COLUMNS)


def write_segment_input(clr, bins, segments, destination, c):
    """One temporary requested-grid, cis-band input; raw observations stay intact."""
    mapping = np.full(len(bins), -1, dtype=np.int64)
    membership = np.full(len(bins), -1, dtype=np.int64)
    frames = []
    cursor = 0
    for group_id, row in enumerate(segments.itertuples(index=False)):
        ids = np.arange(row.native_first_bin, row.native_last_bin + 1)
        frame = bins.iloc[ids][["chrom", "start", "end"]].copy()
        frame["chrom"] = row.segment_id
        frame[["start", "end"]] -= row.start
        frames.append(frame)
        mapping[ids] = np.arange(cursor, cursor + len(ids))
        membership[ids] = group_id
        cursor += len(ids)
    require(frames, "No eligible stripe segments")
    target_bins = pd.concat(frames, ignore_index=True)
    stats = {"raw_contact_sum": 0, "pixel_count": 0}
    starts = bins.start.to_numpy()
    def chunks():
        for frame in pixel_chunks(clr, c["PIXEL_CHUNKSIZE"]):
            require(pd.api.types.is_integer_dtype(frame["count"].dtype), "Stripe caller requires integer raw counts")
            i = frame.bin1_id.to_numpy(); j = frame.bin2_id.to_numpy()
            valid = (mapping[i] >= 0) & (mapping[j] >= 0) & (membership[i] == membership[j])
            valid &= np.abs(starts[j] - starts[i]) <= c["STRIPE_GENOMIC_BELT_BP"]
            count = frame["count"].to_numpy()[valid]
            stats["raw_contact_sum"] += sum(map(int, count))
            stats["pixel_count"] += len(count)
            yield pd.DataFrame(dict(bin1_id=mapping[i[valid]], bin2_id=mapping[j[valid]], count=count))
    cooler.create_cooler(str(destination), target_bins, chunks(), ordered=False, assembly=c["REFERENCE_ID"],
                         dtypes={"count": np.int64})
    require(int(cooler.Cooler(str(destination)).info["sum"]) == stats["raw_contact_sum"], "Stripe temporary raw input violated count conservation")
    return stats


def caller_command(*args):
    return [sys.executable, "-m", "stripepy.main", *map(str, args)]


@contextmanager
def caller_environment():
    """Keep POSIX multiprocessing sockets short without changing task scratch."""
    environment = dict(os.environ, STRIPEPY_NO_TELEMETRY="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    if os.name != "posix":
        yield environment
        return
    # Linux permits 107 address bytes plus NUL; some POSIX systems permit only
    # 103. Reserve 35 bytes for /pymp-XXXXXXXX/listener-XXXXXXXX and the
    # mkdtemp prefix/random suffix before accepting an inherited base directory.
    prefix = "m2c-ipc-"
    parent = None
    for candidate in (*[environment.get(key) for key in ("TMPDIR", "TMP", "TEMP")], "/tmp"):
        if not candidate:
            continue
        path = Path(candidate).resolve()
        if (len(os.fsencode(path)) + 1 + len(prefix) + 8 + 35 <= 103 and
                path.is_dir() and os.access(path, os.W_OK | os.X_OK)):
            parent = path
            break
    require(parent is not None, "StripePy requires a writable short POSIX temporary directory for multiprocessing sockets")
    with tempfile.TemporaryDirectory(prefix=prefix, dir=parent) as temporary:
        os.chmod(temporary, 0o700)
        environment.update(TMPDIR=temporary, TMP=temporary, TEMP=temporary)
        yield environment


def run_caller(command, stdout_path=None):
    with caller_environment() as environment:
        if stdout_path is None:
            subprocess.run(command, check=True, env=environment)
        else:
            with Path(stdout_path).open("w", encoding="utf-8", newline="") as handle:
                subprocess.run(command, check=True, env=environment, stdout=handle)


def canonicalize(native, segment, location, sample, resolution, c):
    """Explicit HDF5 inclusive-bin to BED half-open conversion; no gene strand."""
    left, right, top, bottom, seed = [int(native[k]) for k in ("left_bound", "right_bound", "top_bound", "bottom_bound", "seed")]
    require(location in {"LT", "UT"} and 0 <= left <= seed <= right and 0 <= top <= bottom,
            "Invalid StripePy native geometry")
    require((location == "LT" and top == seed) or (location == "UT" and bottom == seed), "StripePy direction/seed mismatch")
    length = int(segment["end"]) - int(segment["start"])
    require(max(right, bottom) < int(segment["bin_count"]), "StripePy rectangle exceeds its input segment")
    offset = int(segment["start"])
    a0, a1 = offset + left * resolution, offset + min((right + 1) * resolution, length)
    b0, b1 = offset + top * resolution, offset + min((bottom + 1) * resolution, length)
    rc, inner, outer = [float(native[k]) for k in ("rel_change", "inner_mean", "outer_mean")]
    if not np.isfinite([rc, inner, outer]).all() or outer <= 0:
        return None, "undefined_contrast"
    if inner <= outer:
        return None, "nonpositive_contrast"
    if rc < c["STRIPE_MIN_RELATIVE_CHANGE"]:
        return None, "below_relative_change_threshold"
    if bottom == top:
        return None, "zero_excursion"
    direction = "downstream" if location == "LT" else "upstream"
    identity = [c["REFERENCE_ID"], sample["sample_id"], resolution, segment["chrom"], a0, a1, b0, b1, direction]
    row = dict(zip(BEDPE_COLUMNS, [segment["chrom"], a0, a1, segment["chrom"], b0, b1]))
    row.update(stripe_id="stripe_" + digest(identity)[:20], sample_id=sample["sample_id"], condition=sample["condition"],
               reference_id=c["REFERENCE_ID"], resolution=resolution, anchor_chrom=segment["chrom"], anchor_start=a0, anchor_end=a1,
               body_chrom=segment["chrom"], body_start=b0, body_end=b1, seed_bp=offset + seed * resolution,
               direction=direction, matrix_location=location, width_bp=a1-a0, length_bp=b1-b0,
               view=segment["view"], segment_id=segment["segment_id"], native_left_bin=left, native_right_bin=right,
               native_top_bin=top, native_bottom_bin=bottom, relative_change=rc, top_persistence=float(native["top_persistence"]),
               inner_mean=inner, outer_mean=outer, coefficient_of_variation=float(native.get("cfx_of_variation", np.nan)),
               inference_status="descriptive")
    return row, None


def stripes(c, scope):
    from .science import context
    require(c.get("STRIPE_NORMALIZATION", "none") == "none", "Pinned stripe adapter supports raw normalization only")
    require(c.get("STRIPE_MASK_POLICY", "common_segments") == "common_segments", "Stripe adapter supports common_segments mask policy only")
    software = probe()
    samples = validate_metadata(c)["samples"]
    sample = next((s for s in samples if s["sample_id"] == scope), None)
    require(sample is not None, "Stripe scope must be a biological observation")
    root = Path(c["run_dir"]); directory = root / "03_features/stripes"; directory.mkdir(parents=True, exist_ok=True)
    outputs = []
    for resolution in stripe_resolutions(c):
        base = directory / f"{scope}.{resolution}"
        clr, bins, mask, labels, views, _ = context(c, scope, resolution)
        segments = segments_for_context(bins, mask, labels, views, resolution)
        accepted = segments.loc[segments.status == "eligible"]
        # Avoid Path.with_suffix treating the resolution as an extension.
        def product(suffix): return Path(str(base) + suffix)
        segment_path = product(".segments.tsv"); segments.to_csv(segment_path, sep="\t", index=False)
        rows, native_rows, exclusions, commands = [], [], [], []
        stats = {"raw_contact_sum": 0, "pixel_count": 0}
        hdf5 = product(".hdf5"); hdf5.unlink(missing_ok=True)
        product(".caller.tsv").unlink(missing_ok=True)
        if len(accepted):
            scratch = Path(c.get("SCRATCH_ROOT", str(directory))) / "stripe-inputs" / scope
            scratch.mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryDirectory(prefix=f"{resolution}-", dir=scratch) as temporary:
                input_path = Path(temporary) / "segments.cool"
                stats = write_segment_input(clr, bins, accepted, input_path, c)
                command = caller_command("call", input_path, resolution, "--output-file", hdf5,
                    "--genomic-belt", c["STRIPE_GENOMIC_BELT_BP"], "--max-width", c["STRIPE_MAX_WIDTH_BP"],
                    "--nproc", c["THREADS_STRIPES"], "--min-chrom-size", 1, "--low-memory")
                run_caller(command); commands.append(command)
                require(hdf5.is_file(), "StripePy did not produce its HDF5")
            from stripepy.data_structures import ResultFile
            mapping = accepted.set_index("segment_id").to_dict("index")
            with ResultFile(hdf5) as result:
                require(result.resolution == resolution, "StripePy result resolution mismatch")
                for segment_id in result.chromosomes:
                    require(segment_id in mapping, "StripePy returned an unregistered segment")
                    segment = dict(mapping[segment_id], segment_id=segment_id)
                    for location in ("LT", "UT"):
                        # Pinned 1.3.0 cannot read quartiles for an empty triangle.
                        # Geometry is independent of quartiles, so skip only a
                        # demonstrably empty geometry table; other errors fail.
                        if result.get(segment_id, "geo_descriptors", location).empty:
                            continue
                        for native in result.get(segment_id, "stripes", location).to_dict("records"):
                            native_rows.append(dict(segment_id=segment_id, matrix_location=location, **native))
                            row, reason = canonicalize(native, segment, location, sample, resolution, c)
                            if row is not None: rows.append(row)
                            else: exclusions.append(dict(segment_id=segment_id, matrix_location=location, seed_bin=native["seed"], reason=reason))
            pd.DataFrame(native_rows, columns=["segment_id", "matrix_location", "seed", "left_bound", "right_bound", "top_bound", "bottom_bound", "top_persistence", "inner_mean", "outer_mean", "rel_change", "cfx_of_variation"]).to_csv(product(".caller.tsv"), sep="\t", index=False)
            outputs.extend([hdf5, product(".caller.tsv")])
        calls = pd.DataFrame(rows, columns=STRIPE_COLUMNS).drop_duplicates("stripe_id")
        calls.to_csv(product(".tsv"), sep="\t", index=False)
        calls[BEDPE_COLUMNS].to_csv(product(".bedpe"), sep="\t", header=False, index=False)
        # Header-only empty BEDPE is intentionally nonempty for stage validation.
        if calls.empty: product(".bedpe").write_text("# " + "\t".join(BEDPE_COLUMNS) + "\n", encoding="utf-8")
        pd.DataFrame(exclusions, columns=EXCLUSION_COLUMNS).to_csv(product(".exclusions.tsv"), sep="\t", index=False)
        mask_base = root / f"00_metadata/masks/{resolution}"
        atomic_json(product(".json"), dict(software, sample_id=scope, condition=sample["condition"], reference_id=c["REFERENCE_ID"],
            resolution=resolution, input_uri=clr.uri, input_sha256=read_json(root/"00_metadata/inputs.json").get(clr.uri.split("::")[0], {}).get("sha256"),
            mask_sha256=file_hash(str(mask_base)+".npz"), views_sha256=file_hash(str(mask_base)+".views.tsv"),
            commands=commands, raw_input=stats, temporary_input_retained=False, mask_policy="common_segments",
            eligible_segment_count=len(accepted), rejected_segment_count=len(segments)-len(accepted),
            excluded_bin_count=int(((labels >= 0) & ~mask).sum()), candidate_count=len(rows)+len(exclusions), retained_count=len(calls),
            selected_bin_count=int((labels >= 0).sum()), usable_bin_count=int(((labels >= 0) & mask).sum()),
            eligible_bin_count=int(accepted.bin_count.sum()), selected_view_bp=int((views.end-views.start).sum()),
            eligible_segment_bp=int((accepted.end-accepted.start).sum()),
            hdf5_coordinates="segment-local inclusive bin indices; use segments.tsv for original genomic offsets",
            export_coordinates="original-genome 0-based half-open; HDF5 inclusive upper bounds incremented by one bin",
            relative_change_units="native percent: abs(inner_mean-outer_mean)/outer_mean*100",
            positive_contrast_required=True, minimum_relative_change=c["STRIPE_MIN_RELATIVE_CHANGE"],
            native_export="ResultFile geometry-first reader; skips empty triangles to avoid pinned StripePy 1.3.0 empty-quartile assertion",
            min_segment_bins=3, min_caller_chromosome_bp=1, status="detected" if len(calls) else "no_retained_calls",
            absence_interpretation="No calls within eligible contiguous segments does not establish biological absence or genome-wide coverage",
            inference_status="descriptive", biological_p_values="not_computed", telemetry="disabled"))
        outputs.extend([segment_path, product(".tsv"), product(".bedpe"), product(".exclusions.tsv"), product(".json")])
    return outputs


def stripes_summary(c):
    """Dynamic condition inventory and pairwise anchor/direction support, no FDR."""
    root = Path(c["run_dir"]); directory = root/"03_features/stripes"
    samples = validate_metadata(c)["samples"]; summary = []; matches = []
    for res in stripe_resolutions(c):
        tables = {}
        for sample in samples:
            sid = sample["sample_id"]; table = pd.read_csv(directory/f"{sid}.{res}.tsv", sep="\t")
            tables[sid] = table
            summary.append(dict(sample_id=sid, condition=sample["condition"], reference_id=c["REFERENCE_ID"], resolution=res,
                                stripe_count=len(table), biological_replicates_in_condition=sum(s["condition"] == sample["condition"] for s in samples),
                                distinct_subjects_in_condition=len({s.get("subject_id", s["sample_id"]) for s in samples if s["condition"] == sample["condition"]}), inference_status="descriptive"))
        for i, left in enumerate(samples):
            for right in samples[i+1:]:
                indexed = {}
                for key, frame in tables[right["sample_id"]].groupby(["anchor_chrom", "direction"], sort=False):
                    ordered = frame.sort_values("anchor_start", kind="stable")
                    indexed[key] = (ordered, ordered.anchor_start.to_numpy(), int(ordered.width_bp.max()))
                for a in tables[left["sample_id"]].itertuples(index=False):
                    group = indexed.get((a.anchor_chrom, a.direction))
                    if group is None: continue
                    candidates, starts, maximum_width = group
                    lower = np.searchsorted(starts, a.anchor_start-maximum_width, side="right")
                    upper = np.searchsorted(starts, a.anchor_end, side="left")
                    same = candidates.iloc[lower:upper]
                    same = same[same.anchor_end > a.anchor_start]
                    for b in same.itertuples(index=False):
                        overlap = max(0, min(a.body_end, b.body_end)-max(a.body_start, b.body_start))
                        matches.append(dict(resolution=res, sample1=left["sample_id"], sample2=right["sample_id"], stripe1=a.stripe_id, stripe2=b.stripe_id,
                            same_condition=left["condition"] == right["condition"], direction=a.direction, body_overlap_bp=overlap,
                            body_reciprocal_overlap=min(overlap/a.length_bp, overlap/b.length_bp), inference_status="descriptive"))
    sample_path = directory/"summary.tsv"; pd.DataFrame(summary).to_csv(sample_path, sep="\t", index=False)
    match_path = directory/"matches.tsv"
    pd.DataFrame(matches, columns=["resolution","sample1","sample2","stripe1","stripe2","same_condition","direction","body_overlap_bp","body_reciprocal_overlap","inference_status"]).to_csv(match_path, sep="\t", index=False)
    metadata_path = directory/"summary.json"
    atomic_json(metadata_path, dict(reference_id=c["REFERENCE_ID"], resolutions=stripe_resolutions(c),
        condition_members={condition:[s["sample_id"] for s in samples if s["condition"] == condition] for condition in sorted({s["condition"] for s in samples})},
        matching="strict anchor overlap plus same direction, within one resolution; pairwise many-to-many links preserved",
        counts="calls per biological observation; biological_replicates_in_condition counts manifest observations, distinct_subjects_in_condition preserves donor counts; no pooled-map replicates", inference_status="descriptive", biological_p_values="not_computed"))
    return [sample_path, match_path, metadata_path]
