import json
from pathlib import Path

import pytest

from microc2compare.contracts import ContractError, enabled, normalize, read_json
from microc2compare.matrix import resolution_set
from microc2compare.state import prerequisites, relevant_config


def selected(tmp_path, config, names):
    values = {k: str(v).lower() if isinstance(v, bool) else str(v) for k, v in config.items()}
    values.update(M2C_ANALYSES=names, COMPARISONS="contrasts.tsv", RUN_POOL="true", POOLS_MANIFEST="pools.tsv", RUN_DEPTH_SENSITIVITY="true")
    (tmp_path/"contrasts.tsv").write_text("comparison_id\tcase_condition\treference_condition\tdesign\tinference\tmodules\nc_vs_r\tc\tr\t~condition\tdescriptive\tloops;insulation\n")
    return normalize(values, tmp_path)


def test_saddle_requires_only_compartments(tmp_path, config):
    c = selected(tmp_path, config, "saddle")
    assert c["RUN_COMPARTMENTS"] and c["AGGREGATE_MODULES"] == "saddle"
    assert not c["RUN_INSULATION"] and not c["RUN_LOOPS"] and not c["RUN_HICREP"]
    assert resolution_set(c) == [100000]
    assert prerequisites(c, "aggregates") == ["register", "compartments"]
    assert not enabled(c, "quantify")


def test_tads_selection_and_insulation_comparison_do_not_require_loops(tmp_path, config):
    c = selected(tmp_path, config, "tads")
    assert c["RUN_TADS"] and c["RUN_INSULATION"] and not c["RUN_DIFFERENTIAL"]
    assert resolution_set(c) == [10000] and not enabled(c, "quantify")
    c = selected(tmp_path, config, "insulation_comparison")
    assert c["RUN_INSULATION"] and c["RUN_DIFFERENTIAL"]
    assert prerequisites(c, "differential") == ["register", "insulation"]
    assert not enabled(c, "quantify")


def test_optional_derived_analyses_require_explicit_opt_in(tmp_path, config):
    values = {k: str(v).lower() if isinstance(v, bool) else str(v) for k, v in config.items()}
    with pytest.raises(ContractError, match="RUN_POOL=true"):
        normalize(dict(values, M2C_ANALYSES="pool"), tmp_path)
    with pytest.raises(ContractError, match="RUN_DEPTH_SENSITIVITY=true"):
        normalize(dict(values, M2C_ANALYSES="depth_sensitivity"), tmp_path)


def test_threshold_change_preserves_registration_settings(tmp_path, config):
    c = selected(tmp_path, config, "tads")
    changed = dict(c, INSULATION_BOUNDARY_THRESHOLD=.5)
    assert relevant_config(c, "register") == relevant_config(changed, "register")
    assert relevant_config(c, "insulation") != relevant_config(changed, "insulation")
    changed = dict(c, TAD_MIN_LENGTH_BP=100000)
    assert relevant_config(c, "insulation") == relevant_config(changed, "insulation")
    assert relevant_config(c, "tads") != relevant_config(changed, "tads")


def test_invalid_analysis_and_threshold_rejected(tmp_path, config):
    with pytest.raises(ContractError, match="Unknown analysis"):
        selected(tmp_path, config, "guess")
    values = {k: str(v).lower() if isinstance(v, bool) else str(v) for k, v in config.items()}
    with pytest.raises(ContractError, match="finite and nonnegative"):
        normalize(dict(values, TAD_SENSITIVITY_THRESHOLDS="nan"), tmp_path)
