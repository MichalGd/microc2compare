"""Consumer signatures include the parameters and helper code they execute."""
from pathlib import Path

import pytest

from microc2compare.contracts import digest
from microc2compare import state


@pytest.mark.parametrize("key", ["CHROMOSIGHT_MAX_MISSING_FRACTION", "MAX_DENSE_BINS"])
def test_chromosight_discovery_signature_binds_missingness_and_dense_limit(config, key):
    before = state.relevant_config(config, "chromosight_discovery")
    assert before[key] == config[key]
    changed = dict(config)
    changed[key] += .05 if key.endswith("FRACTION") else 1
    after = state.relevant_config(changed, "chromosight_discovery")
    assert after[key] == changed[key]
    assert digest(before) != digest(after)


@pytest.mark.parametrize("module", ["qc_parallel.py", "expected_cache.py", "parallel.py"])
def test_qc_source_signature_changes_when_a_consumed_helper_changes(monkeypatch, module):
    qc_before = state.source_digest("qc")
    unaffected_before = state.source_digest("compartments")
    original = state.file_hash
    touched = []
    def changed_hash(path):
        if Path(path).name == module:
            touched.append(path)
            return "changed-helper-content"
        return original(path)
    monkeypatch.setattr(state, "file_hash", changed_hash)
    assert state.source_digest("qc") != qc_before
    assert touched
    assert state.source_digest("compartments") == unaffected_before
