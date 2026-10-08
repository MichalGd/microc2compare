"""Plot consumer checks use explicit model-result fixtures, not R/calibration evidence."""
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from microc2compare.contracts import ContractError, atomic_json, normalize, read_json
from microc2compare.state import relevant_config
from microc2compare.volcano import prepare_points, render_loop_volcanoes

COMPARISON = dict(comparison_id="treated_vs_control", case_condition="treated", reference_condition="control",
                  inference="formal", modules="loops", design="~condition")


def model_results():
    return pd.DataFrame(dict(feature_id=["higher", "lower", "null", "filtered", "invalid", "zero"],
        logFC=[2., -1., .2, np.nan, 1., 3.], PValue=[.001, .02, .9, np.nan, -.1, 0.],
        FDR=[.01, .04, .95, np.nan, .01, 0.], retained=[True, True, True, False, True, True],
        inference_status=["formal"]*6, validation_status=["uncalibrated"]*6))


def completed_fixture(tmp_path, config, inference="formal", validation="uncalibrated"):
    c = dict(config, run_dir=str(tmp_path))
    comparison = dict(COMPARISON, inference=inference)
    directory = tmp_path/"06_comparisons"/comparison["comparison_id"]
    directory.mkdir(parents=True)
    result = directory/"loops.tsv"
    frame = model_results().assign(inference_status=inference, validation_status=validation)
    frame.to_csv(result, sep="\t", index=False)
    parameters = directory/"parameters.json"
    atomic_json(parameters, dict(inference=inference, validation_status=validation, fdr=.05))
    return c, comparison, result, parameters


def test_model_direction_fdr_family_and_magnitude_display_only():
    original = model_results()
    before = original.copy(deep=True)
    points, counts = prepare_points(original, .05, 1.5)
    assert counts == dict(input_rows=6, not_retained=1, retained_invalid_statistics=1,
        plotted_rows=4, fdr_pass_rows=3, highlighted_rows=2, zero_pvalue_display_clipped=1)
    points = points.set_index("feature_id")
    assert points.loc["higher", "plot_group"] == "Higher in case"
    assert points.loc["lower", "plot_group"] == "FDR passed; not highlighted"
    assert points.loc["null", "plot_group"] == "FDR above threshold"
    assert points.loc["higher", "minus_log10_PValue"] == pytest.approx(3.)
    assert points.loc["zero", "PValue"] == 0 and np.isfinite(points.minus_log10_PValue).all()
    assert points.loc["zero", "pvalue_display_clipped"]
    pd.testing.assert_frame_equal(original, before)
    ordinary, _ = prepare_points(original, .05)
    assert ordinary.set_index("feature_id").loc["lower", "plot_group"] == "Lower in case"
    pd.testing.assert_series_equal(points.FDR, ordinary.set_index("feature_id").FDR)


@pytest.mark.parametrize("column,value", [("logFC", np.inf), ("PValue", 1.1), ("FDR", -1.), ("FDR", np.nan)])
def test_invalid_statistics_are_excluded_with_counts(column, value):
    frame = model_results().iloc[:1].copy()
    frame[column] = value
    points, counts = prepare_points(frame, .05)
    assert points.empty and counts["retained_invalid_statistics"] == 1


def test_missing_and_duplicate_model_identity_fail():
    with pytest.raises(ContractError, match="model columns"):
        prepare_points(model_results().drop(columns="PValue"), .05)
    with pytest.raises(ContractError, match="unique"):
        prepare_points(pd.concat([model_results(), model_results()]), .05)


def test_small_raw_pvalue_does_not_override_model_fdr():
    frame = model_results().iloc[:1].assign(PValue=.001, FDR=.2)
    points, _ = prepare_points(frame, .05)
    assert points.plot_group.iloc[0] == "FDR above threshold"


def test_exports_original_statistics_and_development_label(tmp_path, config, monkeypatch):
    monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path/"plot-cache"))
    c, comparison, result, parameters = completed_fixture(tmp_path, config)
    before = hashlib.sha256(result.read_bytes()).hexdigest()
    outputs, plots, statuses = render_loop_volcanoes(c, [comparison], [result, parameters])
    assert {p.suffix for p in outputs} == {".png", ".svg", ".tsv", ".json"}
    assert all(p.is_file() and p.stat().st_size for p in outputs)
    assert len(plots) == 1 and "UNCALIBRATED DEVELOPMENT" in plots[0][0]
    assert statuses[0]["plotted_rows"] == 4 and statuses[0]["retained_invalid_statistics"] == 1
    points = pd.read_csv(tmp_path/"08_reports/volcano/treated_vs_control.points.tsv", sep="\t")
    assert points.loc[points.feature_id.eq("zero"), "PValue"].iloc[0] == 0
    assert "biological differential PValue" in read_json(tmp_path/"08_reports/volcano/status.json")["axes"]["y"]
    assert hashlib.sha256(result.read_bytes()).hexdigest() == before


@pytest.mark.parametrize("inference", ["formal", "preliminary"])
def test_calibrated_and_preliminary_labels(tmp_path, config, monkeypatch, inference):
    monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path/"plot-cache"))
    c, comparison, result, parameters = completed_fixture(tmp_path, config, inference, "calibrated")
    _, plots, statuses = render_loop_volcanoes(c, [comparison], [result, parameters])
    assert inference.upper() in plots[0][0]
    assert statuses[0]["validation_status"] == "calibrated"


@pytest.mark.parametrize("reason", ["descriptive", "stale", "disabled", "empty"])
def test_no_fake_plot_when_unavailable(tmp_path, config, monkeypatch, reason):
    monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path/"plot-cache"))
    c, comparison, result, parameters = completed_fixture(tmp_path, config)
    current = [result, parameters]
    if reason == "descriptive": comparison["inference"] = "descriptive"
    if reason == "stale": current = []
    if reason == "disabled": c["PLOT_LOOP_VOLCANO"] = False
    if reason == "empty": model_results().iloc[:0].to_csv(result, sep="\t", index=False)
    outputs, plots, statuses = render_loop_volcanoes(c, [comparison], current)
    assert not plots and len(outputs) == 1 and statuses[0]["status"] == "unavailable"
    assert statuses[0]["reason"] and not list((tmp_path/"08_reports/figures").glob("volcano.*"))


def test_status_mismatch_fails_and_unused_comparisons_are_ignored(tmp_path, config):
    c, comparison, result, parameters = completed_fixture(tmp_path, config)
    atomic_json(parameters, dict(inference="formal", validation_status="calibrated", fdr=.05))
    with pytest.raises(ContractError, match="validation_status mismatch"):
        render_loop_volcanoes(c, [comparison], [result, parameters])
    c["_comparison_modules"] = ["insulation"]
    _, plots, statuses = render_loop_volcanoes(c, [comparison], [result, parameters])
    assert not plots and not statuses


def test_display_configuration_is_validated_and_does_not_change_models(tmp_path, config):
    environment = {k: str(v).lower() if isinstance(v, bool) else str(v) for k,v in config.items()}
    for name, value in (("VOLCANO_LABEL_TOP_N", "21"), ("VOLCANO_MIN_ABS_LOG2_EFFECT", "-1")):
        with pytest.raises(ContractError): normalize(dict(environment, **{name: value}), str(tmp_path))
    modified = dict(config, PLOT_LOOP_VOLCANO=False, VOLCANO_MIN_ABS_LOG2_EFFECT=2., VOLCANO_LABEL_TOP_N=0)
    assert relevant_config(config, "differential") == relevant_config(modified, "differential")
    assert relevant_config(config, "quantify") == relevant_config(modified, "quantify")
    assert relevant_config(config, "report") != relevant_config(modified, "report")
    from microc2compare.science import calibration_key
    domain = {"observations": 6, "residual_df": 4}
    assert calibration_key(config, "loops", COMPARISON, {}, domain) == calibration_key(modified, "loops", COMPARISON, {}, domain)


def test_report_uses_committed_comparisons_and_links_exports(tmp_path, config, monkeypatch):
    monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path/"plot-cache"))
    from microc2compare import report
    c, comparison, result, parameters = completed_fixture(tmp_path, config)
    for key in c:
        if key.startswith("RUN_") and isinstance(c[key], bool): c[key] = False
    c["RUN_DIFFERENTIAL"] = c["RUN_FINAL_REPORT"] = True
    manifest = tmp_path/"comparisons.tsv"
    pd.DataFrame([comparison]).to_csv(manifest, sep="\t", index=False)
    c["COMPARISONS"] = str(manifest)
    samples = [dict(sample_id="S1", condition="treated", assay="microc", subject_id="S1")]
    monkeypatch.setattr(report, "validate_metadata", lambda c: dict(samples=samples, comparisons=[comparison]))
    monkeypatch.setattr(report, "complete", lambda c, stage, scope: stage=="differential")
    monkeypatch.setattr(report, "scopes", lambda *a: ["run"])
    atomic_json(tmp_path/"00_metadata/inputs.json", {})
    atomic_json(tmp_path/"completion/run/differential.done.json", dict(stage="differential", scope="run",
        outputs=[dict(path=str(p.relative_to(tmp_path))) for p in (result, parameters)]))
    stale = tmp_path/"06_comparisons/stale/loops.tsv"
    stale.parent.mkdir(); model_results().to_csv(stale, sep="\t", index=False)
    products = report.render(c, "run")
    page = (tmp_path/"08_reports/index.html").read_text(encoding="utf-8")
    assert "UNCALIBRATED DEVELOPMENT" in page and "volcano.treated_vs_control.svg" in page
    assert "treated_vs_control.points.tsv" in page and "stale/loops.tsv" not in page
    assert any(p.suffix==".svg" for p in products)
    # A stale completion cannot authorize plots or link older plot exports.
    monkeypatch.setattr(report, "complete", lambda *a: False)
    report.render(c, "run")
    page = (tmp_path/"08_reports/index.html").read_text(encoding="utf-8")
    assert "stale files are not plotted" in page and "volcano.treated_vs_control.svg" not in page
