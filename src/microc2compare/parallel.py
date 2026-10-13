"""Explicit CPU/RAM bounds for Python workers inside a single stage task."""
from .contracts import allocation_limits, require


def bounded_workers(c, stage, count, worker_memory_mb=None, parent_memory_mb=0):
    from .state import RESERVATIONS
    key = RESERVATIONS[stage]
    limits = allocation_limits()
    threads = min(c["THREADS_" + key], c["TOTAL_CPUS"], max(1, int(limits["effective_cpus"])))
    memory = min(c["MEMORY_MB_" + key], c["MAX_MEMORY_MB"])
    if limits["available_memory_bytes"] is not None:
        memory = min(memory, limits["available_memory_bytes"] // 1000000)
    require(memory > parent_memory_mb, "Stage RAM reservation cannot cover its parent process")
    worker_memory_mb = worker_memory_mb or memory - parent_memory_mb
    require(worker_memory_mb > 0, "Worker RAM reservation must be positive")
    workers = min(count, threads, (memory - parent_memory_mb) // worker_memory_mb,
                  c["MAX_PARALLEL_JOBS"], c["MAX_PARALLEL_SAMPLES"], c["MAX_PARALLEL_MATRIX"])
    if c.get("LIMIT_MATRIX_READERS", False):
        workers = min(workers, c["MAX_PARALLEL_IO"])
    require(workers >= 1, "Stage RAM/CPU budget cannot admit one worker; increase its explicit reservation")
    return int(workers)
