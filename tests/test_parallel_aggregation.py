"""Numerical agreement and recovery on fresh matrices, without old run caches."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

from microc2compare.contracts import atomic_json, file_hash, normalize, read_json
from microc2compare.science import _aggregates_single, aggregates, context, register
from microc2compare.state import initialize
from microc2compare.expected_cache import expected_cached
from microc2compare.quantification import quantify_observations

ROOT = Path(__file__).resolve().parents[1]


def fresh(tmp_path):
    spec = importlib.util.spec_from_file_location('new_aggregate_generator', ROOT/'scripts/generate_demo.py')
    generator = importlib.util.module_from_spec(spec); spec.loader.exec_module(generator)
    data = tmp_path/'fresh'; generator.generate(data, development=True)
    values = {k: str(v['default']).lower() if isinstance(v['default'], bool) else str(v['default'])
              for k, v in read_json(ROOT/'schemas/config.schema.json')['properties'].items()}
    values.update(PROJECT='demo', REFERENCE_ID='demo', SAMPLESHEET='samples.csv', REFERENCE_REGISTRY='references.json',
                  OUTDIR='results', SCRATCH_ROOT='scratch', RUN_LOOPS='false', RUN_INSULATION='false',
                  RUN_COMPARTMENTS='false', RUN_HICREP='false', RUN_DIFFERENTIAL='false', RUN_QC='false',
                  RUN_AGGREGATES='true', AGGREGATE_MODULES='apa', TOTAL_CPUS='4', THREADS_AGGREGATES='4',
                  MAX_PARALLEL_SAMPLES='4', MAX_PARALLEL_MATRIX='4', MAX_PARALLEL_JOBS='4',
                  AGGREGATE_WORKER_MEMORY_MB='1000', QUANTIFY_WORKER_MEMORY_MB='1000', THREADS_QUANTIFY='4')
    c = normalize(values, data); initialize(c, 1, 'fresh-test'); register(c, 'run')
    calls = pd.DataFrame([dict(chrom1='chrDemo', start1=1000000, end1=1010000, chrom2='chrDemo',
                               start2=1200000, end2=1210000, feature_id='candidate_test', representative_resolution=10000)])
    path = Path(c['run_dir'])/'03_features/loops/consensus.tsv'; path.parent.mkdir(parents=True)
    calls.to_csv(path, sep='\t', index=False)
    # Direct pileups use the explicit fixed test universe, independently of discovery.
    c['_apa_candidates'] = calls; c['_apa_resolution'] = 10000
    return c, data, calls


def test_parallel_pileups_match_serial_and_recover_one_corrupted_unit(tmp_path):
    c, data, calls = fresh(tmp_path)
    before = {p.name: file_hash(p) for p in data.glob('*.mcool')}
    _aggregates_single(c, 'run')
    location = Path(c['run_dir'])/'05_aggregates/apa/10000'
    expected = {p.stem: np.load(p)['mean'].copy() for p in location.glob('*.npz')}
    products = aggregates(c, 'run')
    for sid, mean in expected.items():
        np.testing.assert_allclose(np.load(location/(sid+'.npz'))['mean'], mean, rtol=1e-12, atol=1e-12, equal_nan=True)
    times = {sid: (location/(sid+'.npz')).stat().st_mtime_ns for sid in expected}
    c['THREADS_AGGREGATES'] = 2
    aggregates(c, 'run')
    assert times == {sid: (location/(sid+'.npz')).stat().st_mtime_ns for sid in expected}
    (location/'C1.npz').write_bytes(b'corrupted interrupted product')
    aggregates(c, 'run')
    np.testing.assert_allclose(np.load(location/'C1.npz')['mean'], expected['C1'], equal_nan=True)
    assert all((location/(sid+'.npz')).stat().st_mtime_ns == time for sid, time in times.items() if sid != 'C1')
    assert before == {p.name: file_hash(p) for p in data.glob('*.mcool')}


def test_expected_cache_key_includes_mask_and_weights(tmp_path, monkeypatch):
    from microc2compare import matrix
    c, data, calls = fresh(tmp_path)
    args = context(c, 'C1', 10000)
    first = expected_cached(c, 'C1', 10000, *args)
    compute = matrix.expected_table
    counter = []
    def counted(*args): counter.append(1); return compute(*args)
    monkeypatch.setattr(matrix, 'expected_table', counted)
    second = expected_cached(c, 'C1', 10000, *args)
    pd.testing.assert_frame_equal(first, second)
    assert not counter
    clr, bins, mask, labels, views, values = args
    changed = values.copy(); changed[np.flatnonzero(mask)[0]] *= 2
    expected_cached(c, 'C1', 10000, clr, bins, mask, labels, views, changed)
    assert len(counter) == 1
    mask = mask.copy(); mask[np.flatnonzero(mask)[0]] = False
    expected_cached(c, 'C1', 10000, clr, bins, mask, labels, views, changed)
    assert len(counter) == 2


def test_parallel_counts_match_serial_counts_and_exposure(tmp_path):
    from microc2compare.contracts import validate_metadata
    c, data, calls = fresh(tmp_path)
    samples = validate_metadata(c)['samples']
    measurements = {'discovered_loops': (calls, pd.DataFrame())}
    serial = quantify_observations(dict(c, THREADS_QUANTIFY=1), samples, 10000, measurements)
    parallel = quantify_observations(c, samples, 10000, measurements)
    for a, b in zip(serial, parallel):
        av, ae, ar = a['discovered_loops']; bv, be, br = b['discovered_loops']
        pd.testing.assert_series_equal(av, bv); assert ae == be and ar == br


def test_expected_cache_invalidates_registered_core_and_runtime_identity(tmp_path, monkeypatch):
    import copy
    from microc2compare import matrix, state
    c, data, _ = fresh(tmp_path)
    args = context(c, 'C1', 10000)
    environment = Path(c['run_dir'])/'00_metadata/environments.json'
    core = {'identity': 'first-reviewed-core'}
    atomic_json(environment, {'core': core})
    runtime = state.inventory()
    monkeypatch.setattr(state, 'inventory', lambda: copy.deepcopy(runtime))
    first = expected_cached(c, 'C1', 10000, *args)
    record = Path(c['run_dir'])/'00_metadata/expected_cache/10000/C1.json'
    first_identity = read_json(record)['identity_digest']
    compute, recomputed = matrix.expected_table, []
    def counted(*args):
        recomputed.append(1)
        return compute(*args)
    monkeypatch.setattr(matrix, 'expected_table', counted)
    pd.testing.assert_frame_equal(first, expected_cached(c, 'C1', 10000, *args))
    assert not recomputed
    core = {'identity': 'second-reviewed-core'}
    atomic_json(environment, {'core': core})
    pd.testing.assert_frame_equal(first, expected_cached(c, 'C1', 10000, *args))
    second_identity = read_json(record)['identity_digest']
    assert second_identity != first_identity and len(recomputed) == 1
    assert read_json(record)['identity']['registered_core_environment'] == core
    expected_cached(c, 'C1', 10000, *args)
    assert len(recomputed) == 1
    runtime['packages']['numpy'] = 'simulated-test-runtime-version-change'
    pd.testing.assert_frame_equal(first, expected_cached(c, 'C1', 10000, *args))
    assert read_json(record)['identity_digest'] != second_identity and len(recomputed) == 2
    expected_cached(c, 'C1', 10000, *args)
    assert len(recomputed) == 2
