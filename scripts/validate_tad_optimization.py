#!/usr/bin/env python3
"""Replay TADs from a completed pilot into a new directory; compare all products."""
import argparse
import json
import os
import shutil
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
import numpy as np
import pandas as pd
from microc2compare.contracts import atomic_json, file_hash, require
from microc2compare.tads import run


def compare_table(old, new):
    left, right = pd.read_csv(old, sep='\t'), pd.read_csv(new, sep='\t')
    require(list(left.columns) == list(right.columns) and len(left) == len(right), f'Table layout differs: {old}')
    for column in left:
        a, b = left[column], right[column]
        if a.dtype.kind == 'f' and b.dtype.kind == 'f':
            np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-12, equal_nan=True, err_msg=f'{old}: {column}')
        else:
            pd.testing.assert_series_equal(a, b, check_exact=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-run', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--workers', type=int, default=8)
    args = parser.parse_args()
    baseline, output = Path(args.baseline_run).resolve(), Path(args.out).resolve()
    require(not output.exists(), 'Acceptance directory exists; choose a new directory')
    require(args.workers > 0, 'Workers must be positive')
    config_path = baseline/'00_metadata/config.json'
    config = json.loads(config_path.read_text())
    require(Path(config['run_dir']).resolve() == baseline, 'Baseline configuration/run directory differs')
    record = baseline/'completion/run/tads.done.json'
    completed = json.loads(record.read_text())
    require(completed['status'] in ('success', 'success_empty'), 'Baseline TAD stage was not successful')
    inputs = [config_path, record]
    inputs += sorted((baseline/'03_features/insulation').glob('*.tsv'))
    resolution = config['INSULATION_RESOLUTION']
    inputs += [baseline/f'00_metadata/masks/{resolution}.npz', baseline/f'00_metadata/masks/{resolution}.bins.tsv']
    products = []
    for entry in completed['outputs']:
        relative = Path(entry['path'])
        require(not relative.is_absolute() and '..' not in relative.parts, 'Invalid completion path')
        path = baseline/relative
        require(file_hash(path) == entry['sha256'], f'Baseline product hash differs: {path}')
        products.append(path)
    before = {str(path): file_hash(path) for path in inputs+products}
    os.umask(0o077)
    output.mkdir(parents=True)
    for path in inputs[2:]:
        target = output/path.relative_to(baseline)
        target.parent.mkdir(parents=True, exist_ok=True)
        # Only read-only insulation/mask inputs are linked, never output trees.
        try:
            target.symlink_to(path)
        except (OSError, NotImplementedError):
            # Some workstation accounts cannot create symlinks. Copying a TSV
            # or mask is also safe; source matrices are never copied or changed.
            shutil.copyfile(path, target)
    config.update(run_dir=str(output), THREADS_TADS=args.workers,
                  MEMORY_MB_TADS=min(16000, config['MAX_MEMORY_MB']))
    atomic_json(output/'replay-config.json', config)
    started = time.perf_counter()
    generated = run(config, 'run')
    seconds = time.perf_counter()-started
    generated_relative = {path.relative_to(output) for path in generated if path.name != 'performance.json'}
    expected_relative = {path.relative_to(baseline) for path in products if path.name != 'performance.json'}
    require(generated_relative == expected_relative, 'Generated product set differs from the accepted pilot')
    for relative in sorted(expected_relative):
        old, new = baseline/relative, output/relative
        if relative.suffix == '.tsv':
            compare_table(old, new)
        else:
            require(old.read_bytes() == new.read_bytes(), f'Product differs: {relative}')
    for path, expected in before.items():
        require(file_hash(path) == expected, f'Baseline changed during acceptance: {path}')
    old_seconds = None
    benchmark = baseline/'benchmarks/tads.run.json'
    if benchmark.exists():
        old_seconds = json.loads(benchmark.read_text()).get('wall_seconds')
    performance = json.loads((output/'03_features/tads/performance.json').read_text())
    result = dict(status='passed', baseline_run=str(baseline), products_compared=len(expected_relative),
                  float_tolerance=dict(rtol=1e-12, atol=1e-12), coordinates_ids_classes_exclusions='exact',
                  old_wall_seconds=old_seconds, new_replay_wall_seconds=seconds,
                  speedup=old_seconds/seconds if old_seconds else None,
                  workers=performance['workers'], baseline_unchanged=True,
                  timing_note='old stage wall time versus new replay wall time; includes different wrapper overhead')
    atomic_json(output/'acceptance.json', result)
    print(json.dumps(result, indent=2), flush=True)
    print('TAD optimization accepted; completed pilot and matrices were not modified.', flush=True)


if __name__ == '__main__':
    main()
