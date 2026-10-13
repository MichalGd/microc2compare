#!/usr/bin/env python3
"""Bounded synthetic algorithm benchmark; not a whole-genome runtime prediction."""
import argparse
import importlib.util
import json
import platform
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
import numpy as np
import pandas as pd
from microc2compare import features, science, tads
from microc2compare.contracts import atomic_json, file_hash


def timed(function):
    started = time.perf_counter()
    value = function()
    return value, time.perf_counter()-started


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    output = Path(args.out)
    if output.exists():
        raise SystemExit('Refusing to overwrite benchmark results')
    spec = importlib.util.spec_from_file_location('microc2compare._benchmark_reference', ROOT/'tests/fixtures/performance_0_3_2.py')
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    c = {key: prop['default'] for key, prop in json.loads((ROOT/'schemas/config.schema.json').read_text())['properties'].items()}
    c.update(INSULATION_BOUNDARY_THRESHOLD='0.2', LOOP_MIN_DISTANCE_BP=5000, LOOP_MAX_DISTANCE_BP=2000000)
    rng = np.random.default_rng(33)
    n = 19000
    table = pd.DataFrame(dict(chrom=np.repeat([f'chr{i}' for i in range(1, 20)], 1000),
                              region=np.repeat([f'chr{i}' for i in range(1, 20)], 1000),
                              start=np.tile(np.arange(1000)*10000, 19)))
    table['end'] = table.start+10000
    table['boundary_strength_100000'] = np.where(np.arange(n)%11 == 0, rng.uniform(0, 1, n), np.nan)
    table['log2_insulation_score_100000'] = rng.normal(size=n)
    table['valid_pixel_fraction_100000'] = 1.
    mask = np.ones(n, bool)
    a, a_seconds = timed(lambda: old.call_domains(table, mask, c, 100000, 'S', .2, True))
    b, b_seconds = timed(lambda: tads.call_domains(table, mask, c, 100000, 'S', .2, True))
    for left, right in zip(a, b):
        pd.testing.assert_frame_equal(left, right, check_exact=False, atol=1e-12, rtol=1e-12)
    results = {'domain_calling': dict(old_seconds=a_seconds, new_seconds=b_seconds, speedup=a_seconds/b_seconds,
                                     rows=n, domains=len(a[0]), workers=1)}
    index = tads.InsulationIndex(table, mask, c)
    queries = [(str(row.chrom), str(row.region), int(row.start), int(row.end)) for row in a[0].head(800).itertuples()]
    def legacy_evaluable(query):
        chrom, region, start, end = query
        area = table.loc[(table.chrom == chrom) & (table.region == region) & (table.start >= start) & (table.start <= end)]
        outer = area.loc[area.start.isin([start, end])]
        good = np.isfinite(area.log2_insulation_score_100000) & (area.valid_pixel_fraction_100000 >= c['INSULATION_MIN_VALID_PIXEL_FRACTION'])
        return bool(len(area) and len(outer) == 2 and np.isfinite(outer.log2_insulation_score_100000).all() and
                    (outer.valid_pixel_fraction_100000 >= c['INSULATION_MIN_VALID_PIXEL_FRACTION']).all() and
                    mask[area.index].all() and good.mean() >= c['INSULATION_MIN_VALID_PIXEL_FRACTION'])
    a, a_seconds = timed(lambda: [legacy_evaluable(q) for q in queries])
    b, b_seconds = timed(lambda: [index.evaluable(q[:2], q[2], q[3], 100000) for q in queries])
    assert a == b
    results['event_coverage'] = dict(old_seconds=a_seconds, new_seconds=b_seconds, speedup=a_seconds/b_seconds, queries=len(queries), workers=1)
    bins = pd.DataFrame(dict(chrom=['chr1']*2000, start=np.arange(2000)*2000, end=(np.arange(2000)+1)*2000))
    i, j = np.triu_indices(len(bins))
    pixels = pd.DataFrame(dict(bin1_id=i, bin2_id=j, count=rng.integers(1, 25, len(i), dtype=np.int64)))
    choose = rng.choice(np.flatnonzero((j-i >= 3) & (j-i <= 1000)), 800, replace=False)
    candidates = pd.DataFrame(dict(chrom1='chr1', start1=i[choose]*2000, end1=(i[choose]+1)*2000,
                                   chrom2='chr1', start2=j[choose]*2000, end2=(j[choose]+1)*2000,
                                   feature_id=[f'f{k}' for k in range(len(choose))]))
    universes = dict(independent=candidates.iloc[:400].reset_index(drop=True), discovered=candidates.iloc[400:].reset_index(drop=True))
    def chunks(*args):
        for offset in range(0, len(pixels), 250000):
            yield pixels.iloc[offset:offset+250000]
    old.pixel_chunks = features.pixel_chunks = chunks
    mask, labels = np.ones(len(bins), bool), np.zeros(len(bins), int)
    a, a_seconds = timed(lambda: {name: old.legacy_quantify(None, frame, bins, mask, labels, c) for name, frame in universes.items()})
    b, b_seconds = timed(lambda: features.quantify_many(None, universes, bins, mask, labels, c))
    for name in universes:
        pd.testing.assert_series_equal(a[name][0], b[name][0], check_exact=True)
        assert a[name][1:] == b[name][1:]
    results['quantification'] = dict(old_seconds=a_seconds, new_seconds=b_seconds, speedup=a_seconds/b_seconds,
                                     stored_pixels=len(pixels), candidate_sets=2, old_pixel_passes=2, new_pixel_passes=1,
                                     includes_HDF5_IO=False, workers=1)
    ex = pd.DataFrame(dict(region1=np.repeat([f'chr{k}' for k in range(1, 20)], 4000),
                           region2=np.repeat([f'chr{k}' for k in range(1, 20)], 4000), dist=np.tile(np.arange(4000), 19),
                           **{'balanced.avg': np.tile(1/(np.arange(4000)+1), 19)}))
    class Map:
        def extent(self, region): return region[1], region[2]
    def raw(clr, region1, region2, c): return np.ones((region1[2]-region1[1], region2[2]-region2[1]))
    old.bounded_raw = science.bounded_raw = raw
    values = np.ones(4000); labels = np.zeros(4000, int); views = pd.DataFrame(dict(name=['chr1']))
    region1, region2 = ('chr1', 100, 151), ('chr1', 500, 551)
    lookup = science.expected_lookup(ex)
    arguments = (Map(), region1, region2, values, labels, views)
    a, a_seconds = timed(lambda: [old.legacy_oe_window(*arguments, ex, c) for _ in range(200)])
    b, b_seconds = timed(lambda: [science.oe_window(*arguments, lookup, c) for _ in range(200)])
    for left, right in zip(a, b): np.testing.assert_array_equal(left, right)
    results['window_expected_lookup'] = dict(old_seconds=a_seconds, new_seconds=b_seconds, speedup=a_seconds/b_seconds,
                                            windows=200, expected_rows=len(ex), includes_HDF5_IO=False, workers=1)
    atomic_json(output, dict(platform=platform.platform(), python=platform.python_version(),
                             fixture_sha256=file_hash(ROOT/'tests/fixtures/performance_0_3_2.py'), results=results,
                             numerical_agreement='passed', scope='synthetic algorithm-only timings; not target-server wall time'))
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
