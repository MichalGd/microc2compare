"""Report-only volcano plots from completed biological loop comparisons."""
from pathlib import Path

import numpy as np
import pandas as pd

from .contracts import IDENTIFIER, atomic_json, read_json, require


def prepare_points(frame, fdr, minimum_effect=0.):
    """Keep tested rows with usable model statistics; never invent P values."""
    required = {"feature_id", "logFC", "PValue", "FDR", "retained"}
    require(required <= set(frame), f"Loop volcano requires model columns: {sorted(required)}")
    require(frame.feature_id.notna().all() and frame.feature_id.astype(str).str.len().gt(0).all()
            and not frame.feature_id.duplicated().any(), "Loop volcano feature IDs must be nonempty and unique")
    effect = pd.to_numeric(frame.logFC, errors="coerce")
    pvalue = pd.to_numeric(frame.PValue, errors="coerce")
    adjusted = pd.to_numeric(frame.FDR, errors="coerce")
    retained = frame.retained.astype(str).str.lower().eq("true")
    usable = np.isfinite(effect) & np.isfinite(pvalue) & pvalue.between(0, 1) & np.isfinite(adjusted) & adjusted.between(0, 1)
    points = frame.loc[retained & usable].copy()
    points["logFC"] = effect.loc[points.index]
    points["PValue"] = pvalue.loc[points.index]
    points["FDR"] = adjusted.loc[points.index]
    points["pvalue_display_clipped"] = points.PValue.eq(0)
    points["minus_log10_PValue"] = -np.log10(points.PValue.clip(lower=np.finfo(float).tiny))
    points["fdr_pass"] = points.FDR.le(fdr)
    highlighted = points.fdr_pass & points.logFC.abs().ge(minimum_effect) & points.logFC.ne(0)
    points["plot_group"] = "FDR above threshold"
    points.loc[points.fdr_pass, "plot_group"] = "FDR passed; not highlighted"
    points.loc[highlighted & points.logFC.gt(0), "plot_group"] = "Higher in case"
    points.loc[highlighted & points.logFC.lt(0), "plot_group"] = "Lower in case"
    counts = {"input_rows": len(frame), "not_retained": int((~retained).sum()),
              "retained_invalid_statistics": int((retained & ~usable).sum()),
              "plotted_rows": len(points), "fdr_pass_rows": int(points.fdr_pass.sum()),
              "highlighted_rows": int(highlighted.sum()),
              "zero_pvalue_display_clipped": int(points.pvalue_display_clipped.sum())}
    return points, counts


def render_loop_volcanoes(c, comparisons, current_outputs):
    """Consume only configured, currently validated differential output paths."""
    import matplotlib.pyplot as plt

    root = Path(c["run_dir"])
    report = root / "08_reports"
    directory = report / "volcano"
    directory.mkdir(parents=True, exist_ok=True)
    figures = report / "figures"
    figures.mkdir(exist_ok=True)
    current = {Path(p).resolve() for p in current_outputs}
    outputs, plots, statuses = [], [], []
    for comparison in comparisons:
        if "loops" not in comparison["modules"].split(";"): continue
        if "loops" not in c.get("_comparison_modules", ["loops", "insulation"]): continue
        cid = comparison["comparison_id"]
        require(IDENTIFIER.fullmatch(cid), "Unsafe loop comparison ID")
        entry = {"comparison_id": cid, "case": comparison["case_condition"], "reference": comparison["reference_condition"],
                 "status": "unavailable", "reason": "", "inference": comparison["inference"]}
        statuses.append(entry)
        result = root / "06_comparisons" / cid / "loops.tsv"
        parameters = result.parent / "parameters.json"
        if not c.get("PLOT_LOOP_VOLCANO", True):
            entry["reason"] = "Disabled by PLOT_LOOP_VOLCANO=false."
            continue
        if comparison["inference"] == "descriptive":
            entry["reason"] = "Descriptive effects have no biological P values or FDR; no volcano plot is appropriate."
            continue
        if result.resolve() not in current or parameters.resolve() not in current:
            entry["reason"] = "Current validated differential outputs are unavailable; stale files are not plotted."
            continue
        saved = read_json(parameters)
        require(saved["inference"] == comparison["inference"], "Loop volcano inference metadata mismatch")
        require(saved["validation_status"] in ("calibrated", "uncalibrated"), "Loop volcano requires explicit statistical validation status")
        fdr = float(saved["fdr"])
        require(np.isfinite(fdr) and 0 <= fdr <= 1, "Loop volcano FDR threshold is invalid")
        # Preserve literal IDs such as "null" or "NA"; numeric columns are
        # explicitly parsed below, including missing/invalid model values.
        frame = pd.read_csv(result, sep="\t", keep_default_na=False, dtype={"feature_id": str})
        for column, expected in (("inference_status", saved["inference"]), ("validation_status", saved["validation_status"])):
            require(column in frame and frame[column].eq(expected).all(), f"Loop volcano {column} mismatch")
        minimum = c.get("VOLCANO_MIN_ABS_LOG2_EFFECT", 0.)
        points, counts = prepare_points(frame, fdr, minimum)
        entry.update(counts, validation_status=saved["validation_status"], fdr_threshold=fdr,
                     minimum_abs_log2_effect_for_highlighting=minimum,
                     zero_pvalue_display_floor=float(np.finfo(float).tiny),
                     source=str(result.relative_to(root)).replace("\\", "/"))
        if points.empty:
            entry["reason"] = "No retained loops have finite model effects and valid biological P values/FDR."
            continue
        points_path = directory / (cid + ".points.tsv")
        points.to_csv(points_path, sep="\t", index=False)
        outputs.append(points_path)
        fig, ax = plt.subplots(figsize=(7, 5))
        colours = {"FDR above threshold": "#a0a0a0", "FDR passed; not highlighted": "#8a67a8",
                   "Lower in case": "#2878b5", "Higher in case": "#d64d45"}
        for group, colour in colours.items():
            rows = points.loc[points.plot_group.eq(group)]
            if len(rows): ax.scatter(rows.logFC, rows.minus_log10_PValue, s=12, alpha=.65,
                                     c=colour, edgecolors="none", rasterized=True, label=f"{group} (n={len(rows)})")
        ax.axvline(0, color="#777777", lw=.7)
        if minimum > 0:
            for threshold in (-minimum, minimum): ax.axvline(threshold, color="#aaaaaa", lw=.7, ls="--")
        labels = points.loc[points.plot_group.isin(["Lower in case", "Higher in case"])].sort_values(
            ["PValue", "feature_id"], kind="stable").head(c.get("VOLCANO_LABEL_TOP_N", 5))
        for row in labels.itertuples():
            ax.annotate(str(row.feature_id)[:50], (row.logFC, row.minus_log10_PValue),
                        xytext=(4, 5), textcoords="offset points", fontsize=7)
        badge = saved["inference"].upper()
        if saved["validation_status"] == "uncalibrated": badge += " — UNCALIBRATED DEVELOPMENT"
        title = f"{cid}: {badge}"
        ax.set(title=title, xlabel=f"log2 relative contact allocation: {comparison['case_condition']} / {comparison['reference_condition']}",
               ylabel="−log10 biological differential P value")
        ax.legend(fontsize=7, loc="best")
        note = f"Colours use BH FDR ≤ {fdr:g}; highlighting |log2 effect| ≥ {minimum:g}."
        if counts["zero_pvalue_display_clipped"]: note += " Zero P values use the recorded display floor."
        fig.text(.5, .01, note, ha="center", fontsize=8, wrap=True)
        fig.tight_layout(rect=(0, .06, 1, 1))
        for suffix in ("png", "svg"):
            path = figures / f"volcano.{cid}.{suffix}"
            fig.savefig(path, dpi=150)
            outputs.append(path)
            if suffix == "png": plots.append((title, path))
        plt.close(fig)
        entry.update(status="plotted", reason="Biological model P values; caller enrichment statistics are not used.",
                     labelled_rows=len(labels), points=str(points_path.relative_to(report)).replace("\\", "/"))
    status_path = directory / "status.json"
    atomic_json(status_path, {"comparisons": statuses, "axes": {"x": "edgeR logFC (log2 relative allocation, case/reference)",
                "y": "-log10 biological differential PValue"}, "significance": "BH FDR from the complete retained testing family",
                "display_parameters_do_not_change_tests": True})
    outputs.append(status_path)
    return outputs, plots, statuses
