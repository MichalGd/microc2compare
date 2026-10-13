"""Run-owned, identity-checked expected curves; never changes imported maps."""
from pathlib import Path
import hashlib
import os
import tempfile

import numpy as np
import pandas as pd

from .contracts import atomic_json, digest, file_hash, fingerprint, read_json, require


def array_identity(value):
    a = np.ascontiguousarray(value)
    return {"dtype": str(a.dtype), "shape": list(a.shape),
            "sha256": hashlib.sha256(a.tobytes()).hexdigest()}


def expected_identity(c, scope, resolution, clr, mask, labels, views, values):
    from .state import inventory
    source = Path(str(clr.filename).split("::")[0]).resolve()
    registry = Path(c["run_dir"]) / "00_metadata/inputs.json"
    registered = read_json(registry) if registry.is_file() else {}
    previous = next((v for v in registered.values() if Path(v["path"]).resolve() == source), None)
    if previous:
        stat = source.stat()
        require(stat.st_size == previous['size'] and stat.st_mtime_ns == previous['mtime_ns'],
                f"Input changed since registration: {source}")
    # fingerprint verifies size/mtime before reusing a registered SHA256.
    identity = {"matrix": fingerprint(source, previous), "scope": scope,
                "resolution": resolution, "group": str(getattr(clr, "root", "")),
                "mask": array_identity(mask), "labels": array_identity(labels),
                "weights": array_identity(values), "views": views.to_csv(index=False),
                "ignore_diags": c["INSULATION_IGNORE_DIAGS"],
                "implementation": file_hash(Path(__file__).with_name("matrix.py"))}
    environment_file = Path(c['run_dir'])/'00_metadata/environments.json'
    identity['runtime'] = inventory()
    identity['registered_core_environment'] = read_json(environment_file).get('core') if environment_file.is_file() else None
    return digest(identity), identity


def expected_cached(c, scope, resolution, clr, bins, mask, labels, views, values):
    from .matrix import expected_table
    if not c.get("EXPECTED_CACHE", True):
        return expected_table(clr, bins, mask, labels, views, values, c)
    key, identity = expected_identity(c, scope, resolution, clr, mask, labels, views, values)
    directory = Path(c["run_dir"]) / "00_metadata/expected_cache" / str(resolution)
    table = directory / (scope + ".tsv")
    record = directory / (scope + ".json")
    if table.is_file() and record.is_file():
        saved = read_json(record)
        if saved.get("identity_digest") == key and saved.get("table_sha256") == file_hash(table):
            return pd.read_csv(table, sep="\t", float_precision="round_trip", dtype=saved.get('dtypes'))
    result = expected_table(clr, bins, mask, labels, views, values, c)
    directory.mkdir(parents=True, exist_ok=True)
    fd, pending = tempfile.mkstemp(prefix=".expected-", dir=directory)
    os.close(fd)
    try:
        result.to_csv(pending, sep="\t", index=False, na_rep="NA", float_format="%.17g")
        os.replace(pending, table)
    finally:
        Path(pending).unlink(missing_ok=True)
    atomic_json(record, {"identity_digest": key, "identity": identity,
                         "table_sha256": file_hash(table), "normalization": "stored_weights",
                         "dtypes": {name: str(dtype) for name, dtype in result.dtypes.items()}})
    return result
