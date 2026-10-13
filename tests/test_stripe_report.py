from pathlib import Path

import pandas as pd

from microc2compare.contracts import atomic_json
from microc2compare.stripe_report import render


def test_current_stripe_report_keeps_empty_calls_coverage_and_donor_counts(tmp_path):
    c = dict(run_dir=str(tmp_path), RUN_STRIPES=True, STRIPE_RESOLUTIONS="4000")
    directory = tmp_path / "03_features/stripes"
    directory.mkdir(parents=True)
    samples = [dict(sample_id="A", subject_id="shared", culture_id="cultureA", condition="control"),
               dict(sample_id="B", subject_id="shared", culture_id="cultureB", condition="control"),
               dict(sample_id="C", subject_id="other", condition="treated"),
               dict(sample_id="stale", subject_id="stale", condition="treated")]
    current = set()
    for sample in samples:
        sid = sample["sample_id"]
        calls = [] if sid == "C" else [dict(stripe_id="call"+sid, relative_change=150)]
        path = directory/f"{sid}.4000.tsv"
        pd.DataFrame(calls, columns=["stripe_id", "relative_change"]).to_csv(path, sep="\t", index=False)
        provenance = directory/f"{sid}.4000.json"
        atomic_json(provenance, dict(status="detected" if calls else "no_retained_calls", eligible_segment_count=2,
            rejected_segment_count=1, excluded_bin_count=5, candidate_count=len(calls)+2, retained_count=len(calls),
            selected_view_bp=100000, eligible_segment_bp=80000))
        if sid != "stale":
            current.update([path.resolve(), provenance.resolve()])
    matches_path = directory/"matches.tsv"
    pd.DataFrame([dict(resolution=4000, sample1="A", sample2="B", stripe1="callA", stripe2="callB")]).to_csv(matches_path, sep="\t", index=False)
    current.add(matches_path.resolve())
    outputs, plots, sections = render(c, samples, current)
    assert all(Path(p).is_file() and Path(p).stat().st_size for p in outputs)
    assert len(plots) == 1
    inventory = pd.read_csv(tmp_path/"08_reports/stripe_observations.tsv", sep="\t")
    assert inventory.loc[inventory.sample_id.eq("C"), "stripe_count"].iloc[0] == 0
    assert inventory.loc[inventory.sample_id.eq("stale"), "status"].iloc[0] == "current_products_unavailable"
    assert inventory.loc[inventory.sample_id.eq("A"), "eligible_coverage_percent"].iloc[0] == 80
    assert inventory.loc[inventory.sample_id.eq("A"), "distinct_subjects_in_condition"].iloc[0] == 1
    support = pd.read_csv(tmp_path/"08_reports/stripe_support.tsv", sep="\t")
    assert support.supporting_observations_in_condition.tolist() == [2, 2]
    assert support.supporting_subjects_in_condition.tolist() == [1, 1]
    assert not {"PValue", "FDR"} & set(support.columns)
    text = "".join(sections)
    assert "transformed signal" in text and "neither biological absence" in text and "1/1" in text


def test_disabled_stripes_has_no_report_products(tmp_path):
    assert render(dict(run_dir=str(tmp_path), RUN_STRIPES=False), [], set()) == ([], [], [])
    assert not list(tmp_path.iterdir())
