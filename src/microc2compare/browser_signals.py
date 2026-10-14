"""Export already computed compartment, insulation and boundary measurements.

This module does not open matrices or calculate eigenvectors/insulation.  It uses
the ExportContext supplied by the browser-export command; only the Python
standard library is required.
"""
from __future__ import annotations

import hashlib
import math
import re
import statistics


NA = {"", "na", "nan", "none", "null", "<na>"}
AXES = ("chrom", "start", "end")


def _number(value):
    if value is None or str(value).strip().lower() in NA:
        return None
    value = float(value)
    if not math.isfinite(value):
        # Existing workflow browser tracks omit every nonfinite signal.  In
        # particular log2 insulation can be -inf for a zero local signal; keep
        # it in the copied source TSV rather than inventing a finite value.
        return None
    return value


def _boolean(value):
    text = str(value).strip().lower()
    if value is None or text in NA:
        return None
    if text in ("true", "1"):
        return True
    if text in ("false", "0"):
        return False
    raise ValueError(f"Invalid recorded boolean: {value!r}")


def _coord(row):
    return str(row["chrom"]), int(row["start"]), int(row["end"])


def _safe(value):
    value = str(value)
    clean = re.sub(r"[^A-Za-z0-9_.-]", "_", value)
    if clean != value:
        clean += "_" + hashlib.sha256(value.encode()).hexdigest()[:8]
    if not clean or clean in (".", ".."):
        raise ValueError("Empty or unsafe track label")
    return clean


def _close(left, right):
    return left is not None and right is not None and math.isclose(
        left, right, rel_tol=1e-8, abs_tol=1e-10
    )


def _all_or_none(ctx, paths, feature):
    present = [ctx.exists(path) for path in paths]
    if not any(present):
        ctx.add_missing(feature, "No recorded per-observation outputs found")
        return False
    if not all(present):
        raise ValueError(f"Incomplete {feature} observations: " + ", ".join(
            path for path, found in zip(paths, present) if not found
        ))
    return True


def _grid(ctx, resolution):
    path = f"00_metadata/masks/{resolution}.bins.tsv"
    rows = ctx.table(path, required=AXES)
    axes = [_coord(row) for row in rows]
    if not axes or len(axes) != len(set(axes)):
        raise ValueError(f"Empty or duplicate registered coordinates: {path}")
    return axes


def _ordered(rows, axes, source):
    keyed = {_coord(row): row for row in rows}
    if len(keyed) != len(rows) or set(keyed) != set(axes):
        raise ValueError(f"Coordinates differ from the complete registered grid: {source}")
    return [keyed[axis] for axis in axes]


def _mean_rows(axes, values, members):
    """Keep complete-case means separate from available-case summaries."""
    rows = []
    for index, (chrom, start, end) in enumerate(axes):
        present = [(sid, values[sid][index]) for sid in members
                   if values[sid][index] is not None]
        numbers = [value for _, value in present]
        mean = statistics.mean(numbers) if numbers else None
        rows.append({
            "chrom": chrom, "start": start, "end": end,
            "mean_complete": mean if len(numbers) == len(members) else "NA",
            "mean_available": mean if numbers else "NA",
            "sd_available": statistics.stdev(numbers) if len(numbers) > 1 else "NA",
            "n_available": len(numbers), "n_total": len(members),
            "contributing_samples": ";".join(sid for sid, _ in present),
        })
    return rows


def _signal(ctx, path, axes, values, **meta):
    ctx.bedgraph(path, [(*axis, value) for axis, value in zip(axes, values)
                       if value is not None], **meta)


def _bed9(axis, name, color="128,128,128"):
    chrom, start, end = axis
    return chrom, start, end, name, 0, ".", start, end, color


def _ab(value):
    if value is None:
        return "unavailable", "128,128,128"
    if value > 0:
        return "A_like", "213,94,0"
    if value < 0:
        return "B_like", "0,114,178"
    return "zero_unassigned", "128,128,128"


def _mean_exports(ctx, base, axes, values, members, **meta):
    rows = _mean_rows(axes, values, members)
    ctx.write_table(base + ".summary.tsv", rows)
    common = [_number(row["mean_complete"]) for row in rows]
    _signal(ctx, base + ".mean_complete.custom.txt", axes, common,
            **dict(meta, description=meta["description"] +
                   f"; equal-weight mean; all {len(members)} observations required"))
    _signal(ctx, base + ".n_available.custom.txt", axes,
            [row["n_available"] for row in rows],
            **dict(meta, description=meta["description"] +
                   f"; finite observations contributing, maximum {len(members)}"))
    return rows, common


def _export_compartments(ctx):
    sample_ids = [sample["sample_id"] for sample in ctx.samples]
    sources = {sid: f"03_features/compartments/{sid}.tsv" for sid in sample_ids}
    if not _all_or_none(ctx, list(sources.values()), "compartments"):
        return
    resolution = int(ctx.config["COMPARTMENT_RESOLUTION"])
    axes = _grid(ctx, resolution)
    views_path = f"00_metadata/masks/{resolution}.views.tsv"
    views = ctx.table(views_path, required=(*AXES, "name"))
    ctx.copy_source(views_path, "tables/compartments/registered_views.source.tsv")
    if len({row["name"] for row in views}) != len(views):
        raise ValueError("Duplicate compartment view names")
    by_chrom = {}
    for row in views:
        by_chrom.setdefault(row["chrom"], []).append(row)
    memberships = []
    for chrom, start, end in axes:
        matches = [row["name"] for row in by_chrom.get(chrom, [])
                   if int(row["start"]) <= start and end <= int(row["end"])]
        if len(matches) > 1:
            raise ValueError("Overlapping compartment views")
        memberships.append(matches[0] if matches else None)
    scores, component_rows = {}, []
    minimum_correlation = float(ctx.config.get("MIN_PHASING_CORRELATION", 0))
    for sid in sample_ids:
        source = sources[sid]
        rows = _ordered(ctx.table(source, required=(*AXES, "oriented_score", "selected_component")), axes, source)
        phase_source = f"03_features/compartments/{sid}.phasing.json"
        phase = ctx.json(phase_source)
        if not isinstance(phase, list):
            raise ValueError(f"Unexpected phasing QC format: {phase_source}")
        phase_map = {row["region"]: row for row in phase}
        if len(phase_map) != len(phase) or set(phase_map) != {v["name"] for v in views}:
            raise ValueError(f"Phasing QC does not match registered views: {sid}")
        if any(row.get("status") not in ("oriented", "ambiguous") for row in phase):
            raise ValueError(f"Unknown recorded phasing status: {sid}")
        sample_scores, selected_bed = [], []
        for row, axis, region in zip(rows, axes, memberships):
            value = _number(row["oriented_score"])
            selected = row["selected_component"]
            q = phase_map.get(region)
            if q and q["status"] == "oriented":
                component = q["component"]
                corr = _number(q["correlation"])
                if component not in ("E1", "E2", "E3") or corr is None or not minimum_correlation <= abs(corr) <= 1:
                    raise ValueError(f"Invalid accepted phasing metadata: {sid}, {region}")
                if selected != component:
                    raise ValueError(f"Selected vector differs from recorded phasing: {sid}, {region}")
                if component not in row:
                    raise ValueError(f"Selected raw component is missing: {sid}, {component}")
                raw = _number(row[component])
                expected = None if raw is None else raw * (1 if corr >= 0 else -1)
                if (value is None) != (expected is None) or (value is not None and not _close(value, expected)):
                    raise ValueError(f"Oriented score differs from selected vector/sign: {sid}, {axis}")
            else:
                if value is not None or str(selected).strip().lower() not in NA:
                    raise ValueError(f"Finite selected vector in unphased/outside-view bin: {sid}, {axis}")
            sample_scores.append(value)
            if value is not None:
                selected_bed.append(_bed9(axis, str(selected), "90,90,90"))
        scores[sid] = sample_scores
        prefix = f"compartments/individual/{_safe(sid)}"
        meta = dict(analysis="compartments", scope="individual", label=sid,
                    source=source, resolution=resolution)
        _signal(ctx, prefix + ".selected_oriented.custom.txt", axes, sample_scores,
                description=f"{sid}: selected eigenvector oriented to registered phasing reference", **meta)
        for component in ("E1", "E2", "E3"):
            if component in rows[0]:
                _signal(ctx, prefix + f".{component}_raw_unphased.custom.txt", axes,
                        [_number(row[component]) for row in rows],
                        description=f"{sid}: raw {component}, arbitrary sign; not an A/B assignment", **meta)
            else:
                ctx.add_missing(f"{sid} {component}", "Raw component absent from recorded compartment table")
        ctx.bed(prefix + ".AB.custom.txt", [_bed9(axis, *_ab(value))
                for axis, value in zip(axes, sample_scores) if value is not None],
                description=f"{sid}: sign of phased selected vector; A-like positive, B-like negative", **meta)
        ctx.bed(prefix + ".selected_component.custom.txt", selected_bed,
                description=f"{sid}: selected E1/E2/E3 identity; not an additional compartment call", **meta)
        ctx.copy_source(source, f"tables/compartments/{_safe(sid)}.source.tsv")
        ctx.copy_source(phase_source, f"tables/compartments/{_safe(sid)}.phasing.json")
        eigen_source = f"03_features/compartments/{sid}.eigenvalues.tsv"
        if ctx.exists(eigen_source):
            ctx.copy_source(eigen_source, f"tables/compartments/{_safe(sid)}.eigenvalues.tsv")
        else:
            ctx.add_missing(f"{sid} compartment eigenvalues", "Recorded eigenvalue table missing")
        for q in phase:
            component_rows.append(dict(sample_id=sid, **q))
    ctx.write_table("tables/compartments/phasing_summary.tsv", component_rows)
    for condition, members in ctx.conditions.items():
        prefix = f"compartments/conditions/{_safe(condition)}.selected_oriented"
        meta = dict(analysis="compartments", scope="condition", label=condition,
                    source=";".join(sources[sid] for sid in members), resolution=resolution,
                    description=f"{condition}: mean of already oriented selected compartment vectors")
        summary, means = _mean_exports(ctx, prefix, axes, scores, members, **meta)
        for index, row in enumerate(summary):
            present = [scores[sid][index] for sid in members if scores[sid][index] is not None]
            n_a, n_b = sum(v > 0 for v in present), sum(v < 0 for v in present)
            row.update(n_A_like=n_a, n_B_like=n_b, n_zero=sum(v == 0 for v in present),
                       mixed_A_B=bool(n_a and n_b),
                       sign_agreement_fraction=max(n_a, n_b) / len(present) if present else "NA",
                       mean_based_assignment=_ab(means[index])[0])
        # Separate enriched table preserves the main numeric summary without overwrite.
        ctx.write_table(prefix + ".assignment_summary.tsv", summary)
        ctx.bed(prefix + ".mean_based_AB.custom.txt", [
            _bed9(axis, "mean_based_" + _ab(value)[0], _ab(value)[1])
            for axis, value in zip(axes, means) if value is not None
        ], **dict(meta, description=f"{condition}: A/B-like sign of complete-case mean; not consensus membership"))
        agreement = [_number(row["sign_agreement_fraction"]) for row in summary]
        _signal(ctx, prefix + ".sign_agreement.custom.txt", axes, agreement,
                **dict(meta, description=f"{condition}: largest A/B sign fraction among available observations; inspect n_available"))
    ctx.add_missing("Condition mean raw E1/E2/E3",
                    "Raw components have arbitrary signs and can exchange biological identity; only existing selected, phased vectors are averaged. No new component alignment was performed.")
    ctx.notes.append("Compartment condition means use the recorded selected and oriented score, which can select E1, E2 or E3 by view and observation. The vectors retain their original algorithmic scaling; different eigenvalues and selected components can affect amplitude. No additional standardization or rescaling is performed. These descriptive means are not new eigenvectors, pooled compartment calls or quantitative compartment-strength estimates. Raw components are exported per observation with their original signs.")


def _export_insulation(ctx):
    sample_ids = [sample["sample_id"] for sample in ctx.samples]
    sources = {sid: f"03_features/insulation/{sid}.tsv" for sid in sample_ids}
    if not _all_or_none(ctx, list(sources.values()), "insulation and boundaries"):
        return
    resolution = int(ctx.config["INSULATION_RESOLUTION"])
    axes = _grid(ctx, resolution)
    frames = {}
    windows = {int(value) for value in str(ctx.config["INSULATION_WINDOWS_BP"]).split(",") if value.strip()}
    for sid, source in sources.items():
        rows = _ordered(ctx.table(source, required=AXES), axes, source)
        observed = {int(match.group(1)) for key in rows[0]
                    if (match := re.fullmatch(r"log2_insulation_score_(\d+)", key))}
        if not windows.issubset(observed):
            raise ValueError(f"Configured insulation window missing: {sid}")
        windows.update(observed)
        frames[sid] = rows
        ctx.copy_source(source, f"tables/insulation/{_safe(sid)}.source.tsv")
    for sid, rows in frames.items():
        for window in windows:
            for field in (f"log2_insulation_score_{window}", f"boundary_strength_{window}", f"is_boundary_{window}"):
                if field not in rows[0]:
                    raise ValueError(f"Incomplete insulation/boundary window across observations: {sid}, {field}")
    for window in sorted(windows):
        scores, strengths, calls = {}, {}, {}
        for sid, rows in frames.items():
            scores[sid] = [_number(row[f"log2_insulation_score_{window}"]) for row in rows]
            strengths[sid] = [_number(row[f"boundary_strength_{window}"]) for row in rows]
            calls[sid] = [_boolean(row[f"is_boundary_{window}"]) for row in rows]
            for axis, score, strength, called in zip(axes, scores[sid], strengths[sid], calls[sid]):
                if called and (score is None or strength is None):
                    raise ValueError(f"Called boundary lacks finite recorded scores: {sid}, {window}, {axis}")
            meta = dict(analysis="insulation", scope="individual", label=sid,
                        source=sources[sid], resolution=resolution, window=window)
            prefix = f"insulation/window_{window}/individual/{_safe(sid)}"
            _signal(ctx, prefix + ".log2_score.custom.txt", axes, scores[sid],
                    description=f"{sid}: recorded log2 insulation, {window}-bp window", **meta)
            _signal(ctx, prefix + ".boundary_strength.custom.txt", axes, strengths[sid],
                    description=f"{sid}: recorded boundary strength, {window}-bp window", **meta)
            boundary_source = f"03_features/boundaries/{sid}.{window}.tsv"
            expected_calls = {axis: strength for axis, strength, called in zip(axes, strengths[sid], calls[sid]) if called}
            if ctx.exists(boundary_source):
                native = ctx.table(boundary_source, required=(*AXES, f"boundary_strength_{window}"))
                native_calls = {_coord(row): _number(row[f"boundary_strength_{window}"]) for row in native}
                if len(native) != len(native_calls) or set(native_calls) != set(expected_calls) or any(
                        not _close(native_calls[axis], expected_calls[axis]) for axis in native_calls):
                    raise ValueError(f"Called-boundary file disagrees with full insulation table: {boundary_source}")
                ctx.copy_source(boundary_source, f"tables/boundaries/window_{window}/{_safe(sid)}.source.tsv")
            else:
                ctx.add_missing(f"{sid} boundary source {window}", "Native boundary TSV missing; browser intervals use existing is_boundary flags in the full insulation table")
            bedrows, details = [], []
            for axis, strength in expected_calls.items():
                name = f"boundary_{axis[0]}_{axis[1]}_{window}"
                bedrows.append(_bed9(axis, name, "0,114,178"))
                details.append(dict(zip(AXES, axis), feature_id=name, boundary_strength=strength,
                                    feature_id_origin="export_coordinate_label", window_bp=window))
            bprefix = f"boundaries/window_{window}/individual/{_safe(sid)}"
            ctx.bed(bprefix + ".called.custom.txt", bedrows,
                    **dict(meta, analysis="boundaries", description=f"{sid}: recorded called boundary bins, {window}-bp window; BED score is unscored"))
            ctx.write_table(bprefix + ".called.tsv", details,
                            fields=[*AXES, "feature_id", "boundary_strength", "feature_id_origin", "window_bp"])
        for condition, members in ctx.conditions.items():
            source = ";".join(sources[sid] for sid in members)
            meta = dict(analysis="insulation", scope="condition", label=condition,
                        source=source, resolution=resolution, window=window,
                        description=f"{condition}: recorded log2 insulation, {window}-bp window")
            prefix = f"insulation/window_{window}/conditions/{_safe(condition)}"
            _mean_exports(ctx, prefix + ".log2_score", axes, scores, members, **meta)
            _mean_exports(ctx, prefix + ".boundary_strength", axes, strengths, members,
                          **dict(meta, description=f"{condition}: recorded boundary strength at identical bins, {window}-bp window"))
            support_rows, union, unanimous, count_track, eligible_track = [], [], [], [], []
            for index, axis in enumerate(axes):
                called_by = [sid for sid in members if calls[sid][index] is True]
                if not called_by:
                    continue
                available = [sid for sid in members if scores[sid][index] is not None and calls[sid][index] is not None]
                total, n_called, n_available = len(members), len(called_by), len(available)
                name = f"boundary_{axis[0]}_{axis[1]}_w{window}:support={n_called}/{total}"
                bed = _bed9(axis, name, "0,114,178")
                union.append(bed)
                if n_called == total:
                    unanimous.append(bed)
                count_track.append((*axis, n_called))
                eligible_track.append((*axis, n_available))
                support_rows.append(dict(zip(AXES, axis), n_called=n_called, n_available=n_available,
                                         n_total=total, support_fraction_all=n_called / total,
                                         support_fraction_available=n_called / n_available if n_available else "NA",
                                         called_samples=";".join(called_by), available_samples=";".join(available),
                                         window_bp=window, coordinate_rule="exact_native_bin; no clustering"))
            bprefix = f"boundaries/window_{window}/conditions/{_safe(condition)}"
            bmeta = dict(meta, analysis="boundaries")
            ctx.bed(bprefix + ".union.custom.txt", union,
                    **dict(bmeta, description=f"{condition}: union of called boundary bins; exact-coordinate support, no coordinate averaging"))
            ctx.bed(bprefix + f".all_{len(members)}_observations.custom.txt", unanimous,
                    **dict(bmeta, description=f"{condition}: boundary called at identical bin in all {len(members)} observations"))
            ctx.bedgraph(bprefix + ".n_called.custom.txt", count_track,
                         **dict(bmeta, description=f"{condition}: number calling each union boundary, denominator {len(members)} observations"))
            ctx.bedgraph(bprefix + ".n_available.custom.txt", eligible_track,
                         **dict(bmeta, description=f"{condition}: observations with finite insulation and known call status at union boundaries"))
            ctx.write_table(bprefix + ".support.tsv", support_rows,
                            fields=[*AXES, "n_called", "n_available", "n_total", "support_fraction_all", "support_fraction_available",
                                    "called_samples", "available_samples", "window_bp", "coordinate_rule"])
        consensus_source = f"03_features/boundaries/consensus.{window}.tsv"
        if ctx.exists(consensus_source):
            consensus = ctx.table(consensus_source, required=(*AXES, "feature_id", "support"))
            ctx.copy_source(consensus_source, f"tables/boundaries/window_{window}/global_consensus.source.tsv")
            ctx.bed(f"boundaries/window_{window}/global_consensus.custom.txt",
                    [_bed9(_coord(row), f"{row['feature_id']}:support={row['support']}", "204,121,167") for row in consensus],
                    analysis="boundaries", scope="all_samples", label="recorded_global_consensus",
                    source=consensus_source, resolution=resolution, window=window,
                    description=f"Recorded all-observation boundary consensus, {window}-bp window; original IDs and coordinates")
        else:
            ctx.add_missing(f"Global boundary consensus {window}", "Recorded global consensus TSV not found")
    ctx.notes.append("Insulation condition tracks average already recorded log2 scores only at identical coordinates with all observations finite. Native boundary condition tracks use exact-bin call support; neighboring called bins are not clustered or shifted. The separately exported recorded global consensus retains the workflow's original clustering and support.")


def export_signals(ctx):
    """Export available signal modules; partial sample sets fail explicitly."""
    _export_compartments(ctx)
    _export_insulation(ctx)
    ctx.notes.append("Signal tracks contain finite measurements only. NA/NaN and positive/negative infinity remain unchanged in copied source TSVs and contribute no value to browser means; missing measurements are never replaced by zero.")
