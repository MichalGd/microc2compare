"""Bounded raw-count extraction across observations, with one scan per grid."""
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
import os

from .parallel import bounded_workers


def quantify_observation(task):
    c, sid, resolution, measurements = task
    from .science import context
    from .features import quantify_many
    clr, bins, mask, labels, views, values = context(c, sid, resolution)
    result = quantify_many(clr, {name: frame for name, (frame, _) in measurements.items()},
                           bins, mask, labels, dict(c, LOOP_RESOLUTION=resolution))
    print(f"FINISHED raw counts/{resolution}/{sid}", flush=True)
    return result


def quantify_observations(c, samples, resolution, measurements):
    reservation = c.get("QUANTIFY_WORKER_MEMORY_MB", 0) or min(2000, c['MEMORY_MB_QUANTIFY'] - 512)
    workers = bounded_workers(c, "quantify", len(samples), reservation, 512)
    tasks = [(c, s['sample_id'], resolution, measurements) for s in samples]
    print(f"Raw quantification/{resolution}: {workers} worker(s), {len(samples)} observations", flush=True)
    if workers == 1: return [quantify_observation(task) for task in tasks]
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[name] = "1"
    with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        return list(pool.map(quantify_observation, tasks))
