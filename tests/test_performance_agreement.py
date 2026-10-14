"""Compare optimized execution with the frozen, accepted 0.3.2 algorithms."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from microc2compare import features, tads
from microc2compare.science import expected_lookup, oe_window, masked_interval_index, interval_usable


def reference():
    name = 'microc2compare._performance_reference'
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, Path(__file__).parent/'fixtures/performance_0_3_2.py')
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


def scores(seed, size=350):
    rng = np.random.default_rng(seed)
    table = pd.DataFrame(dict(chrom=np.repeat(['chr1', 'chr2'], size),
                              region=np.repeat(['arm1', 'arm2'], size),
                              start=np.tile(np.arange(size)*10000, 2)))
    table['end'] = table.start+10000
    for window in (50000, 100000):
        strength = np.full(len(table), np.nan)
        strength[::9] = rng.uniform(0, 1, len(strength[::9]))
        table[f'boundary_strength_{window}'] = strength
        table[f'is_boundary_{window}'] = np.isfinite(strength) & (strength > .15)
        table[f'log2_insulation_score_{window}'] = rng.normal(size=len(table))
        table[f'valid_pixel_fraction_{window}'] = rng.uniform(.8, 1, len(table))
        table.loc[::47, f'log2_insulation_score_{window}'] = np.nan
        table.loc[::57, f'valid_pixel_fraction_{window}'] = .2
    mask = rng.uniform(size=len(table)) > .01
    return table, mask


@pytest.mark.parametrize('seed', range(4))
@pytest.mark.parametrize('threshold', ['Li', 'Otsu', 0, .2, .5])
def test_domain_calls_against_reference(config, seed, threshold):
    table, mask = scores(seed)
    # Sorting must keep the original mask indices aligned.
    table = table.sample(frac=1, random_state=seed)
    c = dict(config, INSULATION_BOUNDARY_THRESHOLD=threshold)
    for window in (50000, 100000):
        expected = reference().call_domains(table, mask, c, window, 'sample', threshold, True)
        actual = tads.call_domains(table, mask, c, window, 'sample', threshold, True)
        for old, new in zip(expected, actual):
            pd.testing.assert_frame_equal(old, new, check_exact=False, rtol=1e-12, atol=1e-12)


def test_event_coverage_and_nested_overlap_indices(config):
    table, mask = scores(9)
    index = tads.InsulationIndex(table, mask, config)
    for chrom, region in [('chr1', 'arm1'), ('chr2', 'arm2'), ('absent', 'x')]:
        for start, end in [(0, 90000), (10000, 30000), (20000, 75000), (180000, 230000), (3490000, 3500000)]:
            area = table.loc[(table.chrom == chrom) & (table.region == region) & (table.start >= start) & (table.start <= end)]
            outer = area.loc[area.start.isin([start, end])]
            score = 'log2_insulation_score_100000'
            fraction = 'valid_pixel_fraction_100000'
            minimum = config['INSULATION_MIN_VALID_PIXEL_FRACTION']
            good = np.isfinite(area[score]) & (area[fraction] >= minimum)
            old = bool(len(area) and len(outer) == 2 and np.isfinite(outer[score]).all() and
                       (outer[fraction] >= minimum).all() and mask[area.index].all() and good.mean() >= minimum)
            assert index.evaluable((chrom, region), start, end, 100000) == old
    frame = pd.DataFrame(dict(chrom=['chr1']*3, region=['arm1']*3, start=[0, 2, 15], end=[100, 3, 20]))
    overlaps = tads.IntervalIndex([frame])
    for start, end in [(50, 60), (100, 110), (-2, 0), (-2, 1), (3, 15)]:
        assert overlaps.overlaps(('chr1', 'arm1'), start, end) == ((frame.start < end) & (frame.end > start)).any()


def test_full_tad_outputs_serial_parallel_and_reference(tmp_path, config, monkeypatch):
    samples = [dict(sample_id=f'{cond}{i}', condition=cond) for cond in ('A', 'B') for i in range(3)]
    comparison = dict(comparison_id='B_vs_A', reference_condition='A', case_condition='B')
    metadata = dict(samples=samples, comparisons=[comparison])
    old = reference()
    monkeypatch.setattr(old, 'validate_metadata', lambda c: metadata)
    monkeypatch.setattr(tads, 'validate_metadata', lambda c: metadata)
    # Avoid host-specific availability in this small, bounded fixture.
    monkeypatch.setattr(tads, 'allocation_limits', lambda: dict(effective_cpus=4, available_memory_bytes=2_000_000_000))
    roots = []
    for mode in ('reference', 'serial', 'parallel'):
        root = tmp_path/mode
        roots.append(root)
        directory = root/'03_features/insulation'
        directory.mkdir(parents=True)
        table, mask = scores(0, 180)
        for n, sample in enumerate(samples):
            frame, _ = scores(n, 180)
            frame.to_csv(directory/f"{sample['sample_id']}.tsv", sep='\t', index=False)
        meta = root/'00_metadata/masks'
        meta.mkdir(parents=True)
        table[['chrom', 'start', 'end']].to_csv(meta/'10000.bins.tsv', sep='\t', index=False)
        np.savez_compressed(meta/'10000.npz', mask=mask)
        c = dict(config, run_dir=str(root), INSULATION_WINDOWS_BP='50000,100000',
                 INSULATION_BOUNDARY_THRESHOLD='Li', THREADS_TADS=2 if mode == 'parallel' else 1)
        (old.run if mode == 'reference' else tads.run)(c, 'run')
    for path in roots[0].rglob('*'):
        if not path.is_file() or path.suffix not in ('.bed', '.json', '.tsv') or 'insulation' in path.parts or '00_metadata' in path.parts:
            continue
        for root in roots[1:]:
            target = root/path.relative_to(roots[0])
            if path.suffix == '.tsv':
                pd.testing.assert_frame_equal(pd.read_csv(path, sep='\t'), pd.read_csv(target, sep='\t'),
                                              check_exact=False, rtol=1e-12, atol=1e-12)
            else:
                assert path.read_bytes() == target.read_bytes()
    import json
    assert json.loads((roots[2]/'03_features/tads/performance.json').read_text())['workers'] == 2
    # Exercise the deployment replay/comparison gate using a completed baseline.
    from microc2compare.contracts import atomic_json, file_hash
    baseline = roots[0]
    baseline_config = dict(config, run_dir=str(baseline), INSULATION_WINDOWS_BP='50000,100000',
                           INSULATION_BOUNDARY_THRESHOLD='Li')
    atomic_json(baseline/'00_metadata/config.json', baseline_config)
    paths = sorted(p for p in baseline.rglob('*') if p.is_file() and
                   p.suffix in ('.tsv', '.bed', '.json') and 'insulation' not in p.parts and '00_metadata' not in p.parts)
    atomic_json(baseline/'completion/run/tads.done.json', dict(status='success', outputs=[
        dict(path=str(p.relative_to(baseline)), sha256=file_hash(p)) for p in paths]))
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location('replay_gate', root/'scripts/validate_tad_optimization.py')
    replay = importlib.util.module_from_spec(spec); spec.loader.exec_module(replay)
    monkeypatch.setattr(sys, 'argv', ['replay', '--baseline-run', str(baseline), '--out', str(tmp_path/'replay'), '--workers', '2'])
    replay.main()
    acceptance = json.loads((tmp_path/'replay/acceptance.json').read_text())
    assert acceptance['status'] == 'passed' and acceptance['baseline_unchanged'] and acceptance['workers'] == 2


def test_tad_memory_and_allocation_admission(config, monkeypatch):
    table, _ = scores(0)
    tables = {str(i): table for i in range(8)}
    monkeypatch.setattr(tads, 'allocation_limits', lambda: dict(effective_cpus=3.5, available_memory_bytes=None))
    assert tads.worker_count(dict(config, THREADS_TADS=8), tables) == 3
    assert tads.worker_count(dict(config, THREADS_TADS=8, TOTAL_CPUS=2), tables) == 2
    with pytest.raises(Exception, match='memory admission'):
        tads.worker_count(dict(config, MEMORY_MB_TADS=1), tables)


def test_consensus_complete_linkage_and_domain_masks(config):
    rng = np.random.default_rng(6)
    for size in (0, 30, 150):
        frames = {}
        for sample in ('A', 'B', 'pool__P'):
            starts = rng.integers(0, 200, size)*2000
            frame = pd.DataFrame(dict(chrom1=rng.choice(['chr1', 'chr2'], size), start1=starts, end1=starts+2000,
                                      chrom2='chr2', start2=starts+400000, end2=starts+402000))
            frames[sample] = frame
        for tolerance in (0, 1, 3):
            pd.testing.assert_frame_equal(reference().legacy_consensus_loops(frames, 'ref', 2000, tolerance),
                                          features.consensus_loops(frames, 'ref', 2000, tolerance), check_exact=True)
    bins = pd.DataFrame(dict(chrom=['chr1']*100+['chr2']*100, start=np.tile(np.arange(100)*10000, 2)))
    bins['end'] = bins.start+10000
    mask = rng.uniform(size=len(bins)) > .03
    index = masked_interval_index(bins, mask)
    for chrom in ('chr1', 'chr2', 'chrX'):
        for start, end in [(0, 100000), (10000, 30000), (25000, 57000), (1000000, 1100000)]:
            inside = (bins.chrom == chrom) & (bins.start >= start) & (bins.end <= end)
            assert interval_usable(index, chrom, start, end) == bool(inside.any() and mask[inside].all())


def test_raw_quantification_one_pass_exact_integers(config, monkeypatch):
    old = reference()
    rng = np.random.default_rng(5)
    n = 120
    bins = pd.DataFrame(dict(chrom=['chr1']*n, start=np.arange(n)*2000, end=(np.arange(n)+1)*2000))
    i, j = np.triu_indices(n)
    chunk = pd.DataFrame(dict(bin1_id=i, bin2_id=j, count=rng.integers(0, 100, len(i), dtype=np.int64)))
    loops = pd.DataFrame([['chr1', a*2000, (a+1)*2000, 'chr1', b*2000, (b+1)*2000, f'loop{a}_{b}']
                          for a, b in [(0, 0), (1, 10), (3, 40), (5, 21), (6, 60), (0, n)]], columns=features.LOOP_COLUMNS)
    loops.loc[2, 'chrom2'] = 'chrX'
    universes = dict(first=loops, second=loops.iloc[1:4].reset_index(drop=True), empty=loops.iloc[:0].copy())
    mask = np.ones(n, bool); mask[6] = False
    labels = np.zeros(n, int); labels[-10:] = 1
    c = dict(config, LOOP_MIN_DISTANCE_BP=5000, LOOP_MAX_DISTANCE_BP=100000)
    calls = []
    def chunks(*args):
        calls.append(1)
        yield chunk.iloc[:4000]
        yield chunk.iloc[4000:]
    monkeypatch.setattr(old, 'pixel_chunks', chunks)
    monkeypatch.setattr(features, 'pixel_chunks', chunks)
    expected = {name: old.legacy_quantify(None, frame, bins, mask, labels, c) for name, frame in universes.items()}
    calls.clear()
    actual = features.quantify_many(None, universes, bins, mask, labels, c)
    assert len(calls) == 1
    for name in universes:
        pd.testing.assert_series_equal(expected[name][0], actual[name][0], check_exact=True)
        assert expected[name][1:] == actual[name][1:]


def test_expected_lookup_oe_exact_zero_nan_distances(config, monkeypatch):
    old = reference()
    class Matrix:
        def extent(self, region):
            return region[1], region[2]
    clr = Matrix()
    values = np.array([1., 2., np.nan, 1., .5, 2.])
    labels = np.array([0, 0, 0, 1, 1, -1])
    views = pd.DataFrame(dict(name=['left', 'right']))
    expected = pd.DataFrame(dict(region1=['left']*3+['right']*2, region2=['left']*3+['right']*2,
                                 dist=[0, 1, 4, 0, 1], **{'balanced.avg': [0., 2., 4., np.nan, 3.]}))
    def raw(clr, a, b, c):
        return np.arange(36).reshape(6, 6)[a[1]:a[2], b[1]:b[2]]
    monkeypatch.setattr(old, 'bounded_raw', raw)
    monkeypatch.setattr('microc2compare.science.bounded_raw', raw)
    lookup = expected_lookup(expected)
    for a, b in [(('chr1', 0, 6), ('chr1', 0, 6)), (('chr1', 0, 3), ('chr1', 1, 5)), (('chr1', 3, 5), ('chr1', 3, 5))]:
        reference_array = old.legacy_oe_window(clr, a, b, values, labels, views, expected, config)
        np.testing.assert_array_equal(oe_window(clr, a, b, values, labels, views, lookup, config), reference_array)
    assert expected_lookup(lookup) is lookup
