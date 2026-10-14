"""Audit-local caching preserves exact identities and fresh corruption checks."""
from collections import Counter
from contextlib import nullcontext
from contextvars import copy_context
import json
from pathlib import Path
import time

import pytest

from microc2compare import state
from microc2compare.contracts import atomic_json, file_hash, fingerprint, read_json


@pytest.fixture
def dependency_graph(tmp_path, config, monkeypatch):
    c = dict(config, run_dir=str(tmp_path / "run"), REFERENCE_ID="fresh_cache_reference")
    for key in c:
        if key.startswith("RUN_"):
            c[key] = False
    graph = {
        "preflight": [], "register": ["preflight"], "qc": ["register"],
        "loops": ["register", "qc"],
        "consensus": ["register", "qc", "loops"],
        "quantify": ["register", "loops", "consensus"],
        "report": ["register", "qc", "loops", "consensus", "quantify"],
    }
    members = [{"sample_id": f"Fresh{i + 1}"} for i in range(3)]
    metadata_calls = []
    def metadata(current):
        metadata_calls.append(current)
        return {"samples": list(members)}
    monkeypatch.setattr(state, "prerequisites", lambda current, stage: list(graph[stage]))
    monkeypatch.setattr(state, "enabled", lambda current, stage: stage in graph)
    monkeypatch.setattr(state, "validate_metadata", metadata)
    root = Path(c["run_dir"])
    input_path = tmp_path / "immutable-input.txt"
    input_path.write_text("independent input identity\n")
    atomic_json(root / "00_metadata/inputs.json", {str(input_path): fingerprint(input_path)})
    atomic_json(root / "00_metadata/environments.json", {"core": {"identity": "synthetic-runtime"}})
    for stage in graph:
        for scope in state.scopes(c, stage):
            product = root / "products" / f"{stage}.{scope}.tsv"
            product.parent.mkdir(parents=True, exist_ok=True)
            product.write_text(f"stage\tscope\tvalue\n{stage}\t{scope}\t0.125\n")
            state.publish_completion(c, stage, scope, [product])
    metadata_calls.clear()
    return c, members, input_path, metadata_calls


def test_exact_signature_and_complete_match_uncached_dependency_audit(dependency_graph, monkeypatch):
    c, _, _, _ = dependency_graph
    cached = state.signature(c, "report", "run")
    assert state.complete(c, stage="report", scope="run")
    with monkeypatch.context() as uncached:
        uncached.setattr(state, "_validation_context", lambda current: nullcontext(None))
        assert state.signature(c, "report", "run") == cached
        assert state.complete(c, "report", "run")
    assert state._VALIDATION_AUDIT.get() is None


@pytest.mark.parametrize("damage", ["output", "dependency", "input", "source", "legacy_record"])
def test_next_top_level_audit_detects_damage(dependency_graph, monkeypatch, damage):
    c, _, input_path, _ = dependency_graph
    assert state.complete(c, "report", "run")
    root = Path(c["run_dir"])
    if damage == "output":
        product = root / "products/qc.Fresh1.tsv"
        original = product.read_bytes()
        product.write_bytes(original.replace(b"0.125", b"9.125"))  # unchanged size, changed hash
    elif damage == "dependency":
        state.record_path(c, "qc", "Fresh1").unlink()
    elif damage == "input":
        input_path.write_text(input_path.read_text() + "modified\n")
    elif damage == "source":
        real_hash = state.file_hash
        monkeypatch.setattr(state, "file_hash", lambda path: "changed-source" if Path(path).name == "scientific.py" else real_hash(path))
    else:
        path = state.record_path(c, "qc", "Fresh1")
        record = read_json(path)
        record.pop("signature")
        atomic_json(path, record)
    assert not state.complete(c, "report", "run")
    assert state._VALIDATION_AUDIT.get() is None


def test_source_digest_standalone_is_fresh_after_an_audit(dependency_graph, monkeypatch):
    c, _, _, _ = dependency_graph
    first = state.source_digest("report")
    assert state.complete(c, "report", "run")
    real_hash = state.file_hash
    monkeypatch.setattr(state, "file_hash", lambda path: "new-source-bytes" if Path(path).name == "scientific.py" else real_hash(path))
    assert state.source_digest("report") != first
    assert state._VALIDATION_AUDIT.get() is None


def test_different_config_nested_audit_does_not_borrow_cached_results(dependency_graph, monkeypatch):
    c, _, _, _ = dependency_graph
    changed = dict(c, REFERENCE_ID="another_reference")
    original = state.relevant_config
    nested = []
    def relevant(current, stage):
        if current is c and stage == "report":
            parent = state._VALIDATION_AUDIT.get()
            nested.append(state.complete(changed, "register", "run"))
            assert state._VALIDATION_AUDIT.get() is parent
        return original(current, stage)
    monkeypatch.setattr(state, "relevant_config", relevant)
    assert state.complete(c, "report", "run")
    assert nested == [False]
    assert not state.complete(changed, "report", "run")
    assert state._VALIDATION_AUDIT.get() is None


def test_scopes_are_defensive_and_metadata_refreshes_next_audit(dependency_graph):
    c, members, _, calls = dependency_graph
    with state._validation_context(c):
        first = state.scopes(c, "qc")
        first.append("caller_mutation")
        assert state.scopes(c, "qc") == [m["sample_id"] for m in members]
        assert state.scopes(c, "loops") == [m["sample_id"] for m in members]
    assert len(calls) == 1
    assert state.complete(c, "report", "run")
    members.append({"sample_id": "NewUnregisteredObservation"})
    assert not state.complete(c, "report", "run")


def test_dependency_graph_hash_work_is_bounded_per_audit(dependency_graph, monkeypatch):
    c, _, _, metadata_calls = dependency_graph
    real_hash = state.file_hash
    calls = Counter()
    def counted(path):
        calls[str(path)] += 1
        return real_hash(path)
    monkeypatch.setattr(state, "file_hash", counted)
    start = time.perf_counter()
    assert state.complete(c, "report", "run")
    cached_seconds = time.perf_counter() - start
    cached = dict(calls)
    assert len(metadata_calls) == 1
    # Source paths are hashed once across all stages. Every unique required
    # product is still hashed during this top-level audit.
    outputs = set()
    root = Path(c["run_dir"])
    for path in (root / "completion").glob("*/*.done.json"):
        outputs.update(str(root / row["path"]) for row in read_json(path)["outputs"])
    assert outputs <= set(cached)
    assert all(count == 1 for count in cached.values())
    calls.clear()
    metadata_calls.clear()
    with monkeypatch.context() as uncached:
        uncached.setattr(state, "_validation_context", lambda current: nullcontext(None))
        start = time.perf_counter()
        assert state.complete(c, "report", "run")
        uncached_seconds = time.perf_counter() - start
    assert sum(calls.values()) > sum(cached.values()) * 10
    assert len(metadata_calls) > 10
    print(json.dumps({"cached_hash_calls": sum(cached.values()), "uncached_hash_calls": sum(calls.values()),
                      "cached_seconds": cached_seconds, "uncached_seconds": uncached_seconds}))


def test_failed_missing_completion_is_not_retained_after_publication(dependency_graph):
    c, _, _, _ = dependency_graph
    record = state.record_path(c, "report", "run")
    record.unlink()
    assert not state.complete(c, "report", "run")
    product = Path(c["run_dir"]) / "products/report.run.tsv"
    state.publish_completion(c, "report", "run", [product])
    assert state.complete(c, "report", "run")


def test_copied_context_cannot_reuse_an_ended_audit(dependency_graph):
    c, _, _, _ = dependency_graph
    with state._validation_context(c):
        assert state.complete(c, "report", "run")
        inherited = copy_context()
    product = Path(c["run_dir"]) / "products/qc.Fresh1.tsv"
    product.write_text("changed after the audit\n")
    assert not inherited.run(state.complete, c, "report", "run")
    assert state._VALIDATION_AUDIT.get() is None
