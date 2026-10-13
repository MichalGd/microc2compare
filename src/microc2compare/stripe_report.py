"""Render descriptive stripe coverage and support from current saved products."""
from pathlib import Path

import numpy as np
import pandas as pd

from .contracts import atomic_json, read_json


def render(c, samples, current):
    """Return report outputs, plot descriptors and HTML sections for fine_report."""
    if not c.get("RUN_STRIPES"):
        return [], [], []
    root = Path(c["run_dir"])
    directory = root / "03_features/stripes"
    report = root / "08_reports"
    figures = report / "figures"
    report.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    current = {Path(p).resolve() for p in current}
    sections = ["<h2>Stripes: coverage and descriptive support</h2><p>StripePy 1.3.0 calls raw contact patterns within contiguous common usable segments. Inherited balancing eligibility defines the mask, but balancing weights do not scale the caller's contacts. Stripes crossing masked gaps are excluded. No retained calls establishes neither biological absence nor complete genome coverage. Each grid is a separate sensitivity measurement.</p>",
                "<p>The narrow anchor and extended body retain upstream/downstream genomic excursion independently of gene strand. Relative change is the caller's percentage abs(inner_mean - outer_mean) / outer_mean × 100 on transformed signal; retained calls also require inner_mean &gt; outer_mean. These percentages are neither contact-count fold changes nor biological differential P values/FDR. Anchor overlap and matching direction provide descriptive geometric support; potentially many-to-many matches do not establish independent replication or regulation.</p>"]
    outputs, plots, inventory, support = [], [], [], []
    by_id = {s["sample_id"]: s for s in samples}
    condition_members = {condition: [s for s in samples if s["condition"] == condition]
                         for condition in sorted({s["condition"] for s in samples})}
    match_path = directory / "matches.tsv"
    matches = pd.read_csv(match_path, sep="\t") if match_path.resolve() in current else None
    for resolution in sorted(int(r) for r in c["STRIPE_RESOLUTIONS"].split(",")):
        for sample in samples:
            sid = sample["sample_id"]
            path, provenance = directory/f"{sid}.{resolution}.tsv", directory/f"{sid}.{resolution}.json"
            if not {path.resolve(), provenance.resolve()} <= current:
                inventory.append(dict(sample_id=sid, condition=sample["condition"], resolution=resolution,
                                      status="current_products_unavailable", inference_status="descriptive"))
                continue
            calls = pd.read_csv(path, sep="\t")
            info = read_json(provenance)
            members = condition_members[sample["condition"]]
            selected_bp, eligible_bp = info.get("selected_view_bp"), info.get("eligible_segment_bp")
            values = pd.to_numeric(calls.relative_change, errors="coerce")
            inventory.append(dict(sample_id=sid, subject_id=sample.get("subject_id", sid),
                experimental_unit_id=sample.get("experimental_unit_id") or sample.get("culture_id", ""),
                condition=sample["condition"], resolution=resolution, stripe_count=len(calls), status=info["status"],
                observations_in_condition=len(members), distinct_subjects_in_condition=len({s.get("subject_id", s["sample_id"]) for s in members}),
                eligible_segments=info["eligible_segment_count"], too_short_segments=info["rejected_segment_count"],
                selected_view_bp=selected_bp, eligible_segment_bp=eligible_bp,
                eligible_coverage_percent=100*eligible_bp/selected_bp if selected_bp else np.nan,
                masked_bins=info["excluded_bin_count"], rejected_predictions=info["candidate_count"]-info["retained_count"],
                relative_change_min_percent=float(values.min()) if len(values) else np.nan,
                relative_change_max_percent=float(values.max()) if len(values) else np.nan, inference_status="descriptive"))
            if matches is None:
                continue
            pairs = matches.loc[matches.resolution.eq(resolution)]
            for call in calls.itertuples(index=False):
                participating = set()
                for pair in pairs.loc[pairs.stripe1.eq(call.stripe_id) | pairs.stripe2.eq(call.stripe_id)].itertuples(index=False):
                    other = pair.sample2 if pair.sample1 == sid else pair.sample1
                    if other in by_id and by_id[other]["condition"] == sample["condition"]:
                        participating.add(other)
                participants = {sid} | participating
                support.append(dict(stripe_id=call.stripe_id, sample_id=sid, condition=sample["condition"], resolution=resolution,
                    supporting_observations_in_condition=len(participants), observations_in_condition=len(members),
                    supporting_subjects_in_condition=len({by_id[p].get("subject_id", p) for p in participants}),
                    distinct_subjects_in_condition=len({s.get("subject_id", s["sample_id"]) for s in members}), inference_status="descriptive"))
    frame = pd.DataFrame(inventory)
    inventory_path = report / "stripe_observations.tsv"
    frame.to_csv(inventory_path, sep="\t", index=False)
    outputs.append(inventory_path)
    sections.append("<p>The table includes every registered observation and configured stripe grid, including empty calls and unavailable current products. Coverage is the eligible segment length divided by selected registered-view length; short usable segments and masked bins remain explicit. Observation denominators and distinct subject counts are shown separately.</p>" + frame.to_html(index=False, na_rep="unavailable"))
    conditions = []
    if "stripe_count" in frame:
        for (resolution, condition), group in frame.groupby(["resolution", "condition"], sort=True):
            counts = group.stripe_count.dropna()
            conditions.append(dict(condition=condition, resolution=resolution,
                registered_observations=len(group), observations_with_current_results=len(counts),
                observations_with_calls=int((counts > 0).sum()), calls_sum=int(counts.sum()),
                mean_calls_per_current_observation=float(counts.mean()) if len(counts) else np.nan,
                inference_status="descriptive"))
    condition_path = report / "stripe_conditions.tsv"
    condition_frame = pd.DataFrame(conditions, columns=["condition", "resolution", "registered_observations", "observations_with_current_results", "observations_with_calls", "calls_sum", "mean_calls_per_current_observation", "inference_status"])
    condition_frame.to_csv(condition_path, sep="\t", index=False)
    outputs.append(condition_path)
    if len(condition_frame):
        sections.append("<p>Condition totals summarize saved calls. The mean uses current observation counts, including zeros. Detection depends on depth and usable coverage; these summaries are descriptive and perform no differential stripe test.</p>" + condition_frame.to_html(index=False, na_rep="unavailable"))
    if matches is not None:
        support_path = report / "stripe_support.tsv"
        support_frame = pd.DataFrame(support, columns=["stripe_id", "sample_id", "condition", "resolution", "supporting_observations_in_condition", "observations_in_condition", "supporting_subjects_in_condition", "distinct_subjects_in_condition", "inference_status"])
        support_frame.to_csv(support_path, sep="\t", index=False)
        outputs.append(support_path)
        sections.append("<p>The downloadable stripe support table counts the focal observation plus other same-condition observations whose calls overlap its anchor and share direction at this grid. It preserves distinct subjects when multiple cultures share a donor. A singleton has support 1/1 without replicate confirmation. Direct pairwise matches are not collapsed into a transitive consensus.</p>")
    else:
        sections.append("<p>Current pairwise matching products are unavailable; replicate support is not calculated from old results.</p>")
    if "stripe_count" in frame:
        import matplotlib.pyplot as plt
        palette = plt.get_cmap("tab20")
        colors = {condition: palette(index % 20) for index, condition in enumerate(condition_members)}
        for resolution, group in frame.groupby("resolution", sort=True):
            group = group.loc[group.stripe_count.notna()]
            if group.empty:
                continue
            fig, ax = plt.subplots(figsize=(9, max(3, min(20, .35*len(group)+1))))
            ax.barh(np.arange(len(group)), group.stripe_count, color=[colors[k] for k in group.condition])
            ax.set_yticks(np.arange(len(group)), labels=[f"{s} ({k})" for s, k in zip(group.sample_id, group.condition)])
            ax.invert_yaxis()
            ax.set(xlabel="Retained stripe calls per observation", title=f"Descriptive stripe inventory ({resolution} bp)")
            fig.tight_layout()
            for suffix in ("png", "svg"):
                picture = figures / f"stripes.{resolution}.counts.{suffix}"
                fig.savefig(picture, dpi=150, bbox_inches="tight")
                outputs.append(picture)
                if suffix == "png":
                    plots.append((f"Stripe calls at {resolution} bp", picture))
            plt.close(fig)
    metadata_path = report / "stripe_report.json"
    atomic_json(metadata_path, dict(inference_status="descriptive", biological_p_values="not_computed",
        coverage="eligible contiguous segment bp / selected registered view bp; no whole-genome absence inference",
        support="focal plus directly anchor-overlapping same-direction calls in same condition and resolution; distinct subjects retained",
        relative_change="native transformed-signal abs(inner_mean-outer_mean)/outer_mean*100; positive contrast required"))
    outputs.append(metadata_path)
    return outputs, plots, sections
