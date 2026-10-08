"""Descriptive comparison outputs from fixed quantified-table fixtures."""
import json

import numpy as np
import pandas as pd
import pytest

from microc2compare import science
from microc2compare.contracts import atomic_json, file_hash
from microc2compare.features import LOOP_COLUMNS


@pytest.mark.parametrize("empty_loops", [False, True])
def test_descriptive_modules_write_parameters_and_preserve_missing_values(tmp_path, config, monkeypatch, empty_loops):
    ids = ["C1", "C2", "C3", "T1", "T2", "T3"]
    samples = [{"sample_id": sid, "condition": "control" if sid.startswith("C") else "treated"} for sid in ids]
    comparison = {"comparison_id": "treated_vs_control", "case_condition": "treated",
                  "reference_condition": "control", "design": "~condition",
                  "inference": "descriptive", "modules": "loops;insulation"}
    monkeypatch.setattr(science, "validate_metadata", lambda c: {"samples": samples, "comparisons": [comparison]})
    def forbidden(*args, **kwargs):
        raise AssertionError("Descriptive comparisons must not launch R")
    monkeypatch.setattr(science, "run_r", forbidden)
    c = dict(config, run_dir=str(tmp_path), REFERENCE_ID="demo")
    atomic_json(tmp_path / "00_metadata/environments.json", {})
    counts = pd.DataFrame([["chrDemo", i * 10000, (i + 1) * 10000,
                            "chrDemo", 100000, 110000, "f" + str(i)] for i in range(3)], columns=LOOP_COLUMNS)
    for sid, count in zip(ids, [10, 20, 30, 40, 80, 120]):
        counts[sid] = [count, 0, np.nan]
    if empty_loops:
        counts = counts.iloc[:0].copy()
    exposure = pd.DataFrame({"sample_id": ids, "eligible_cis_exposure": [100, 200, 300, 200, 400, 600]})
    inputs = []
    for kind in ("discovered_loops", "independent_loops"):
        inputs.append(science.save(counts, tmp_path / f"04_quantification/{kind}/counts.tsv"))
        inputs.append(science.save(exposure, tmp_path / f"04_quantification/{kind}/exposures.tsv"))
    window = c["INSULATION_PRIMARY_WINDOW_BP"]
    for sid in ids:
        case = sid.startswith("T")
        frame = pd.DataFrame({"chrom": ["chrDemo"] * 3, "start": [0, 10000, 20000],
                              "end": [10000, 20000, 30000],
                              f"log2_insulation_score_{window}": [1 if case else 0, 1 if case else np.nan, 0],
                              f"is_boundary_{window}": [True, False, False]})
        inputs.append(science.save(frame, tmp_path / f"03_features/insulation/{sid}.tsv"))
    before = {str(path): file_hash(path) for path in inputs}
    outputs = science.differential(c, "run")
    assert len(outputs) == 6 and all(path.is_file() for path in outputs)
    directory = tmp_path / "06_comparisons/treated_vs_control"
    parameters = json.loads((directory / "parameters.json").read_text())
    assert parameters["design"]["residual_df"] == 4
    assert parameters["design"]["contrast"] == [0, 1]
    assert parameters["validation_status"] == "descriptive"
    loops = pd.read_csv(directory / "loops.tsv", sep="\t")
    if empty_loops:
        assert loops.empty and "log2_allocation_ratio" in loops.columns
    else:
        assert loops.loc[0, "log2_allocation_ratio"] == pytest.approx(1.0)
        assert loops.loc[0, "case_mean_allocation"] == pytest.approx(0.2)
        assert loops.loc[0, "reference_mean_allocation"] == pytest.approx(0.1)
        assert loops.loc[1:, "log2_allocation_ratio"].isna().all()
        assert loops.loc[1, "case_mean_allocation"] == 0
        assert pd.isna(loops.loc[2, "case_mean_allocation"])
    insulation = pd.read_csv(directory / "insulation.tsv", sep="\t")
    np.testing.assert_allclose(insulation.case_minus_reference_log2_insulation, [1, np.nan, 0], equal_nan=True)
    for table in (loops, insulation):
        assert table.PValue.isna().all() and table.FDR.isna().all()
        assert table.inference_status.eq("descriptive").all()
        assert table.validation_status.eq("descriptive").all()
    assert before == {str(path): file_hash(path) for path in inputs}
    assert not list(tmp_path.rglob("*.cool")) and not list(tmp_path.rglob("*.mcool"))
