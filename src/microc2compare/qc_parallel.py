"""Per-grid QC workers with independent matrix handles and bounded admission."""
from concurrent.futures import ProcessPoolExecutor
import math
import multiprocessing
import os
from pathlib import Path

from .contracts import atomic_json, require, validate_metadata
from .parallel import bounded_workers


def resolutions(c):
    from .resolutions import chromosight_discovery_resolutions, loop_resolutions, quantification_resolutions
    grids = {c["PS_RESOLUTION"]}
    if c["RUN_LOOPS"]:
        grids.update(loop_resolutions(c))
    if c.get("RUN_CHROMOSIGHT_DISCOVERY"):
        grids.update(chromosight_discovery_resolutions(c))
    if c["RUN_CHROMOSIGHT"]:
        grids.update(quantification_resolutions(c))
    return sorted(grids)


def _grid_info(c, scope, resolution):
    import cooler
    from .matrix import matrix_uri
    from .resolutions import sample_matrix
    sample = next(s for s in validate_metadata(c)["samples"] if s["sample_id"] == scope)
    return cooler.Cooler(matrix_uri(sample_matrix(c, sample, resolution), resolution)).info


def worker_memory_mb(c, scope, grids):
    """Conservative reservation estimate, not an operating-system RSS limit."""
    override = c.get("QC_WORKER_MEMORY_MB", 0)
    require(isinstance(override, int) and override >= 0, "QC_WORKER_MEMORY_MB must be a nonnegative integer")
    if override:
        return override
    estimates = []
    for resolution in grids:
        info = _grid_info(c, scope, resolution)
        actual_chunk = min(c["PIXEL_CHUNKSIZE"], info["nnz"])
        estimates.append(math.ceil(256 + int(info["nbins"])/1000 + actual_chunk*128/1000000))
    return max(estimates)


def _native_threads_one():
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[name] = "1"


def _grid_job(payload):
    c, scope, resolution = payload
    from .science import context, save
    from .expected_cache import expected_cached
    from .matrix import raw_summary
    clr, bins, mask, labels, views, values = context(c, scope, resolution)
    summary = raw_summary(clr, bins, mask, labels, c)
    contacts = Path(c["run_dir"])/f"01_qc/contacts/{scope}.{resolution}.json"
    atomic_json(contacts, dict(summary, sample_id=scope, resolution=resolution, normalization="inherited_upstream_weights"))
    expected = expected_cached(c, scope, resolution, clr, bins, mask, labels, views, values)
    curve = save(expected, Path(c["run_dir"])/f"01_qc/ps/{scope}.{resolution}.tsv")
    return [contacts, curve]


def run(c, scope):
    if not c["RUN_QC"]:
        path = Path(c["run_dir"])/f"01_qc/hicrep/inputs/{scope}.json"
        atomic_json(path, dict(sample_id=scope, status="raw_only_HiCRep_input_validated",
                              resolution=c["HICREP_RESOLUTION"], normalization="raw_counts_no_balancing"))
        return [path]
    grids = resolutions(c)
    count = bounded_workers(c, "qc", len(grids), worker_memory_mb(c, scope, grids), parent_memory_mb=512)
    jobs = [(c, scope, resolution) for resolution in grids]
    if count == 1:
        return [output for job in jobs for output in _grid_job(job)]
    # Spawn imports the scientific libraries before its initializer runs. Set
    # limits before spawning so imported native pools also receive the budget.
    names = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")
    previous = {name: os.environ.get(name) for name in names}
    _native_threads_one()
    try:
        with ProcessPoolExecutor(max_workers=count, mp_context=multiprocessing.get_context("spawn"),
                                 initializer=_native_threads_one) as pool:
            return [output for result in pool.map(_grid_job, jobs) for output in result]
    finally:
        for name, value in previous.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value
