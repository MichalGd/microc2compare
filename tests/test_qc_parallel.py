"""Actual independent grid workers preserve raw summaries and expected curves."""
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd
import pytest

from microc2compare.contracts import atomic_json, ContractError, file_hash, normalize, read_json, SAMPLE_FIELDS, write_table
from microc2compare import qc_parallel
from microc2compare.science import qc, register
from microc2compare.state import initialize


def project(tmp_path):
    path = tmp_path/"fresh.mcool"
    balancing = {}
    for resolution in (4000, 8000):
        bins = cooler.binnify(pd.Series({"newChr": 320000}), resolution)
        i, j = np.triu_indices(len(bins))
        counts = ((i+1)*3 + j*5).astype(np.int64)
        cooler.create_cooler(f"{path}::/resolutions/{resolution}", bins,
            pd.DataFrame(dict(bin1_id=i, bin2_id=j, count=counts)), ordered=True, mode="a", assembly="newAssembly")
        with h5py.File(path, "r+") as handle:
            values = np.linspace(1, 1.5, len(bins)); values[80000//resolution] = np.nan
            handle[f"resolutions/{resolution}/bins"].create_dataset("weight", data=values)
        balancing[str(resolution)] = dict(weight_column="weight", convention="multiplicative", converged=True,
                                          recipe={"kind": "synthetic_declared_weights"})
    provenance = tmp_path/"provenance.json"
    atomic_json(provenance, dict(reference_id="newAssembly", count_units="unique_contacts",
        filters=dict(mapq=30, min_distance_bp=0, chromosome_policy="newChr"), balancing=balancing,
        dedup_scope="joint_library", library_ids=["L"], source_ids=["source"]))
    write_table(tmp_path/"samples.tsv", [dict(sample_id="fresh", project="newProject", condition="oneCondition",
        assay="microc", biological_replicate="one", subject_id="donor", reference_id="newAssembly",
        matrix_path=str(path), provenance_path=str(provenance))], SAMPLE_FIELDS)
    membership = dict(sample_id="fresh", library_id="L", source_id="source", source_kind="independent_library",
        dedup_group_id="G", dedup_scope="joint_library", source_checksum=file_hash(path))
    write_table(tmp_path/"membership.tsv", [membership], list(membership))
    (tmp_path/"sizes.tsv").write_text("newChr\t320000\n")
    (tmp_path/"views.bed").write_text("newChr\t0\t320000\tchosen\n")
    atomic_json(tmp_path/"reference.json", dict(reference_id="newAssembly", chromsizes="sizes.tsv",
        chromsizes_sha256=file_hash(tmp_path/"sizes.tsv"), views="views.bed", views_sha256=file_hash(tmp_path/"views.bed"),
        chromosome_classes={"autosomes": ["newChr"]}))
    c = normalize(dict(PROJECT="newProject", REFERENCE_ID="newAssembly", SAMPLESHEET="samples.tsv",
        TECHNICAL_MEMBERSHIP="membership.tsv", REFERENCE_REGISTRY="reference.json", OUTDIR="results", SCRATCH_ROOT="scratch",
        M2C_ANALYSES="qc,loops", LOOP_CALLER="mustache", LOOP_RESOLUTION="4000", LOOP_DISCOVERY_RESOLUTIONS="4000,8000",
        LOOP_QUANTIFICATION_RESOLUTIONS="4000,8000", LOOP_MIN_DISTANCE_BP="4000", LOOP_MAX_DISTANCE_BP="320000",
        PS_RESOLUTION="8000", APA_FLANK_BP="80000", TOTAL_CPUS="2", MAX_MEMORY_MB="2000",
        THREADS_QC="2", MEMORY_MB_QC="1800", PIXEL_CHUNKSIZE="100", ENVIRONMENT_MODE="development"), tmp_path)
    initialize(c, 1, "qc-test")
    register(c, "run")
    return c, path


def test_actual_serial_parallel_grids_cache_and_input_identity(tmp_path):
    c, source = project(tmp_path)
    before = file_hash(source)
    assert qc_parallel.resolutions(c) == [4000, 8000]
    c["EXPECTED_CACHE"] = False
    c["THREADS_QC"] = 1
    serial = qc(c, "fresh")
    summaries = {p.name: read_json(p) for p in serial if p.suffix == ".json"}
    curves = {p.name: pd.read_csv(p, sep="\t", float_precision="round_trip") for p in serial if p.suffix == ".tsv"}
    for path in serial:
        path.unlink()
    c["THREADS_QC"] = 2
    assert qc_parallel.bounded_workers(c, "qc", 2, qc_parallel.worker_memory_mb(c, "fresh", [4000, 8000]), parent_memory_mb=512) == 2
    parallel = qc(c, "fresh")
    assert [p.name for p in serial] == [p.name for p in parallel]
    for path in parallel:
        if path.suffix == ".json":
            assert summaries[path.name] == read_json(path)
        else:
            pd.testing.assert_frame_equal(curves[path.name], pd.read_csv(path, sep="\t", float_precision="round_trip"))
    c["EXPECTED_CACHE"] = True
    qc(c, "fresh")
    cache = Path(c["run_dir"])/"00_metadata/expected_cache"
    saved = {p: (p.stat().st_mtime_ns, file_hash(p)) for p in cache.rglob("*") if p.is_file()}
    qc(c, "fresh")
    assert saved == {p: (p.stat().st_mtime_ns, file_hash(p)) for p in cache.rglob("*") if p.is_file()}
    for path in parallel:
        if path.suffix == ".tsv":
            pd.testing.assert_frame_equal(curves[path.name], pd.read_csv(path, sep="\t", float_precision="round_trip"))
    assert before == file_hash(source)
    assert not list(Path(c["run_dir"]).rglob("*.cool")) and not list(Path(c["run_dir"]).rglob("*.mcool"))


def test_hicrep_only_qc_preserves_raw_input_contract(tmp_path):
    c = dict(run_dir=str(tmp_path), RUN_QC=False, HICREP_RESOLUTION=10000)
    outputs = qc(c, "arbitrary")
    assert len(outputs) == 1
    assert read_json(outputs[0]) == dict(sample_id="arbitrary", status="raw_only_HiCRep_input_validated",
                                        resolution=10000, normalization="raw_counts_no_balancing")


def test_native_axis_and_actual_chunk_admission_refuses_before_calculation(tmp_path, monkeypatch):
    c, _ = project(tmp_path)
    c["PIXEL_CHUNKSIZE"] = 100000000
    monkeypatch.setattr(qc_parallel, "_grid_info", lambda *args: {"nbins": 7500000, "nnz": 2000000000})
    assert qc_parallel.worker_memory_mb(c, "fresh", [4000, 8000]) == 20556
    def forbidden(*args):
        raise AssertionError("Insufficient reservation admitted a QC calculation")
    monkeypatch.setattr(qc_parallel, "_grid_job", forbidden)
    with pytest.raises(ContractError, match="cannot admit one worker"):
        qc(c, "fresh")
    c["QC_WORKER_MEMORY_MB"] = 9000
    assert qc_parallel.worker_memory_mb(c, "fresh", [4000, 8000]) == 9000
    with pytest.raises(ContractError, match="cannot admit one worker"):
        qc(c, "fresh")


def test_tiny_actual_headers_fit_one_gb_without_capping_required_memory(tmp_path):
    c, _ = project(tmp_path)
    c["MEMORY_MB_QC"] = 1000
    c["PIXEL_CHUNKSIZE"] = 100000000
    assert qc_parallel.worker_memory_mb(c, "fresh", [4000, 8000]) == 257
    assert qc_parallel.bounded_workers(c, "qc", 2, 257, parent_memory_mb=512) == 1
    assert len(qc(c, "fresh")) == 4
