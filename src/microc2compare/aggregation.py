"""Parallel per-observation pileups with atomic, independently verified units."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import hashlib
import multiprocessing
import os
import tempfile

import numpy as np
import pandas as pd

from .contracts import atomic_json, digest, file_hash, read_json, require, validate_metadata
from .parallel import bounded_workers


def atomic_npz(path, **arrays):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, pending = tempfile.mkstemp(prefix=".pileup-", suffix=".npz", dir=path.parent)
    os.close(fd)
    try:
        np.savez_compressed(pending, **arrays)
        os.replace(pending, path)
    finally:
        Path(pending).unlink(missing_ok=True)


def _unit_key(c, sid, module, resolution):
    from .state import relevant_config
    from .science import context
    from .expected_cache import expected_identity
    clr, bins, mask, labels, views, values = context(c, sid, resolution)
    expected_key, _ = expected_identity(c, sid, resolution, clr, mask, labels, views, values)
    root = Path(c["run_dir"])
    features = {}
    if module == "apa":
        frame = c.get("_apa_candidates")
        if frame is not None:
            features["mapped_candidates"] = hashlib.sha256(frame.to_csv(index=False).encode()).hexdigest()
        for name in ("03_features/loops/consensus.tsv", "03_features/loops/catalogue.tsv"):
            if (root / name).is_file(): features[name] = file_hash(root / name)
    if module == "ada":
        name = "03_features/boundary_derived_domains/consensus.tsv"
        if (root / name).is_file(): features[name] = file_hash(root / name)
    if module == "saddle":
        features["partition"] = hashlib.sha256(c["_saddle_partition"][0].to_csv(index=False).encode()).hexdigest()
    config = relevant_config(c, "aggregates")
    # Include measurement geometry omitted by legacy stage signatures.
    config.update(ADA_RESOLUTION=c["ADA_RESOLUTION"], APA_FLANK_BP=c["APA_FLANK_BP"])
    implementation = {p: file_hash(Path(__file__).with_name(p)) for p in
                      ("aggregation.py", "science.py", "matrix.py", "expected_cache.py")}
    return digest({"sample": sid, "module": module, "resolution": resolution,
                   "expected": expected_key, "features": features, "config": config,
                   "implementation": implementation})


def _unit(c, sid, module, resolution):
    from .science import _aggregates_single
    root = Path(c["run_dir"])
    directory = root / "00_metadata/checkpoints/aggregates" / module / str(resolution) / sid
    config = dict(c, _aggregate_sample_id=sid, _aggregate_unit_dir=str(directory))
    key = _unit_key(config, sid, module, resolution)
    checkpoint = directory / "completed.json"
    if c.get("AGGREGATE_CHECKPOINTS", True) and checkpoint.is_file():
        record = read_json(checkpoint)
        paths = [root / r["path"] for r in record.get("outputs", [])]
        if record.get("identity_digest") == key and paths and all(
                p.is_file() and p.stat().st_size == r["size"] and file_hash(p) == r["sha256"]
                for p, r in zip(paths, record["outputs"])):
            print(f"REUSED {module}/{resolution}/{sid}", flush=True)
            return paths
    checkpoint.unlink(missing_ok=True)
    products = _aggregates_single(config, sid)
    entries = []
    for p in products:
        p = Path(p)
        require(p.is_file() and p.stat().st_size > 0, f"Aggregate unit did not publish {p}")
        entries.append({"path": str(p.relative_to(root)), "size": p.stat().st_size, "sha256": file_hash(p)})
    atomic_json(checkpoint, {"identity_digest": key, "outputs": entries, "sample_id": sid,
                             "module": module, "resolution": resolution, "status": "success"})
    print(f"FINISHED {module}/{resolution}/{sid}", flush=True)
    return [Path(p) for p in products]


def run_units(c, module, resolution):
    """Merge only committed unit tables; each worker opens its own Cooler handles."""
    from .science import save, saddle_partition
    root = Path(c["run_dir"])
    samples = validate_metadata(c)["samples"]
    config = dict(c)
    shared_outputs = []
    if module == "saddle":
        partition, labels = saddle_partition(c, samples)
        require((partition.partition >= 0).any(), "Saddle has no eligible partition")
        config["_saddle_partition"] = (partition, labels)
        shared_outputs.append(save(partition, root / "05_aggregates/saddle/partition.tsv"))
        info = root / "05_aggregates/saddle/partition.json"
        atomic_json(info, {"source": c["SADDLE_PARTITION_SOURCE"], "groups": 50,
                          "partition_hash": digest(partition[["chrom", "start", "end", "partition"]].to_dict("records")),
                          "inference_status": "descriptive"})
        shared_outputs.append(info)
    reservation = c.get("AGGREGATE_WORKER_MEMORY_MB", 0) or min(4000, c['MEMORY_MB_AGGREGATES'] - 512)
    workers = bounded_workers(c, "aggregates", len(samples), reservation, 512)
    print(f"Aggregation {module}/{resolution}: {workers} worker(s), {len(samples)} observations", flush=True)
    products = {}
    if workers == 1:
        for sample in samples:
            sid = sample["sample_id"]
            products[sid] = _unit(config, sid, module, resolution)
    else:
        # spawn avoids inherited HDF5 handles and BLAS thread pools on Linux as well.
        for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
            os.environ[name] = "1"
        with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn")) as pool:
            pending = {pool.submit(_unit, config, s["sample_id"], module, resolution): s["sample_id"] for s in samples}
            try:
                for future in as_completed(pending): products[pending[future]] = future.result()
            except BaseException:
                for future in pending: future.cancel()
                raise
    outputs = list(shared_outputs)
    tables = []
    for sample in samples:
        paths = products[sample["sample_id"]]
        outputs += [p for p in paths if p.suffix != ".tsv"]
        tables += [pd.read_csv(p, sep="\t") for p in paths if p.suffix == ".tsv"]
    merged = pd.concat(tables, ignore_index=True)
    target = root / (f"05_aggregates/apa/{resolution}/summary.tsv" if module == "apa" else
                     "04_quantification/domain_scores/strength.tsv" if module == "ada" else
                     "05_aggregates/saddle/summary.tsv")
    outputs.append(save(merged, target))
    return outputs
