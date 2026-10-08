import json
import os
import sys
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd
import pytest

from microc2compare.contracts import atomic_json,file_hash,normalize,read_json,validate_metadata,write_table
from microc2compare.science import (balance_new_map,cleanup_derived,depth_sensitivity,inspect_inputs,
                                   pools,register,saddle_partition,validate_pool_sum)
from microc2compare.state import (complete,initialize,publish_completion,record_path,restore_required_pools)

REPO=Path(__file__).resolve().parents[1]


@pytest.fixture
def prepared(tmp_path):
    import importlib.util
    spec=importlib.util.spec_from_file_location('generator',REPO/'scripts/generate_demo.py')
    generator=importlib.util.module_from_spec(spec);spec.loader.exec_module(generator)
    root=tmp_path/'data';generator.generate(root,development=True,raw_only=True)
    schema=read_json(REPO/'schemas/config.schema.json')
    values={k:str(v['default']).lower() if isinstance(v['default'],bool) else str(v['default']) for k,v in schema['properties'].items()}
    values.update(PROJECT='demo',ANALYSIS_SET_ID='synthetic_demo',RUN_ID='run001',REFERENCE_ID='demo',ENVIRONMENT_MODE='development',PYTHON_BOOTSTRAP=sys.executable,
                  REFERENCE_REGISTRY='references.json',SAMPLESHEET='samples.csv',OUTDIR='results',SCRATCH_ROOT='scratch',TOTAL_CPUS='2',MAX_MEMORY_MB='2000',
                  RUN_COMPARTMENTS='false',RUN_INSULATION='false',RUN_LOOPS='false',RUN_AGGREGATES='false',RUN_DIFFERENTIAL='false')
    c=normalize(values,str(root));initialize(c,os.getpid(),'test')
    publish_completion(c,'preflight','run',inspect_inputs(c,writing=True))
    publish_completion(c,'register','run',register(c,'run'))
    return root,c


def test_resume_rejects_corruption_and_source_changes(prepared,monkeypatch):
    _,c=prepared
    assert complete(c,'register','run')
    bins=Path(c['run_dir'])/'00_metadata/masks/10000.bins.tsv'
    bins.write_text(bins.read_text()+'corrupt\n')
    assert not complete(c,'register','run')
    publish_completion(c,'register','run',register(c,'run'))
    assert complete(c,'register','run')
    monkeypatch.setattr('microc2compare.state.source_digest',lambda s:'different_code')
    assert not complete(c,'register','run')


def test_mixed_tads_and_selected_insulation_ignore_unused_formal_loop_design(prepared, monkeypatch):
    from types import SimpleNamespace
    root, original = prepared
    (root/'comparisons.tsv').write_text(
        'comparison_id\tcase_condition\treference_condition\tdesign\tinference\tmodules\n'
        'unused_loops\ttreated\tcontrol\t~condition+missing_covariate\tformal\tloops\n'
        'chosen_insulation\ttreated\tcontrol\t~condition\tdescriptive\tinsulation\n')
    c = dict(original, RUN_TADS=True, RUN_INSULATION=True, RUN_DIFFERENTIAL=True, RUN_HICREP=False,
             _comparison_modules=['insulation'])
    # Preflight inspects actual original matrices; no feature calculation runs.
    monkeypatch.setitem(sys.modules, 'cooltools', SimpleNamespace(eigs_cis=None, insulation=None, dots=None))
    requested = []
    def environment(c, kind):
        requested.append(kind)
        assert kind == 'core', 'An unused formal loop comparison requested R'
        return {'test_preflight_environment': kind}
    monkeypatch.setattr('microc2compare.science.environment_record', environment)
    metadata = validate_metadata(c)
    assert len(metadata['comparisons']) == 2  # both contrasts remain available to descriptive TADs
    inspected = inspect_inputs(c)
    assert inspected and requested == ['core']


def test_parameter_changes_reuse_registration_and_invalidate_consumers(prepared):
    from microc2compare.science import qc
    _, original = prepared
    c = dict(original, RUN_INSULATION=True, RUN_TADS=True)
    root = Path(c['run_dir'])
    admission = root/'00_metadata/state-test-admission.json'
    atomic_json(admission, {'purpose': 'Completion dependency test; not scientific insulation validation'})
    publish_completion(c, 'preflight', 'run', [admission])
    assert complete(c, 'register', 'run')  # same imported axes and required grid
    publish_completion(c, 'qc', 'C1', qc(c, 'C1'))
    for sample in validate_metadata(c)['samples']:
        product = root/f"03_features/insulation/{sample['sample_id']}.dependency-test.tsv"
        product.parent.mkdir(parents=True, exist_ok=True)
        product.write_text('purpose\ncompletion dependency fixture\n')
        publish_completion(c, 'insulation', sample['sample_id'], [product])
    tad_product = root/'03_features/tads/dependency-test.tsv'
    tad_product.parent.mkdir(parents=True, exist_ok=True)
    tad_product.write_text('purpose\ncompletion dependency fixture\n')
    publish_completion(c, 'tads', 'run', [tad_product])
    assert complete(c, 'tads', 'run')
    changed = dict(c, TAD_MIN_LENGTH_BP=c['TAD_MIN_LENGTH_BP']+10000)
    publish_completion(changed, 'preflight', 'run', [admission])
    assert complete(changed, 'register', 'run') and complete(changed, 'qc', 'C1')
    assert all(complete(changed, 'insulation', s['sample_id']) for s in validate_metadata(c)['samples'])
    assert not complete(changed, 'tads', 'run')
    changed = dict(c, INSULATION_BOUNDARY_THRESHOLD='0.5')
    publish_completion(changed, 'preflight', 'run', [admission])
    assert complete(changed, 'register', 'run') and complete(changed, 'qc', 'C1')
    assert not complete(changed, 'insulation', 'C1') and not complete(changed, 'tads', 'run')


def test_exact_imported_pool_reuse_without_construction(prepared,monkeypatch):
    root,c=prepared;samples=validate_metadata(c)['samples'];members=samples[:3]
    path=root/'upstream_pool.cool'
    cooler.merge_coolers(str(path),[s['matrix_path']+'::/resolutions/10000' for s in members],mergebuf=100000)
    with h5py.File(path,'r+') as handle:handle['bins'].create_dataset('weight',data=np.ones(len(handle['bins/start'])))
    prov=read_json(members[0]['provenance_path'])
    prov.update(member_sample_ids=[s['sample_id'] for s in members],member_hashes={s['sample_id']:file_hash(s['matrix_path']) for s in members},summation='raw_integer_sum')
    provenance=root/'pool.json';atomic_json(provenance,prov)
    row={'pool_id':'control','pool_scope':'condition','condition':'control','member_sample_ids':'C1;C2;C3','reference_id':'demo','matrix_path':str(path),'provenance_path':str(provenance)}
    manifest=root/'pools.tsv';write_table(manifest,[row],list(row))
    c=dict(c,RUN_POOL=True,POOL_MODE='import',POOLS_MANIFEST=str(manifest))
    before=file_hash(path)
    def forbidden(*a,**kw):raise AssertionError('Imported pool attempted creation/balance')
    monkeypatch.setattr(cooler,'merge_coolers',forbidden);monkeypatch.setattr(cooler,'balance_cooler',forbidden)
    products=pools(c,'run');assert len(products)==1 and file_hash(path)==before
    assert not list((Path(c['run_dir'])/'02_pools').rglob('*.cool'))
    # Equal total alone cannot validate a pool: mutate two pixels oppositely.
    with h5py.File(path,'r+') as handle:
        counts=handle['pixels/count'];counts[0]=counts[0]+1;counts[1]=counts[1]-1
    with pytest.raises(ValueError,match='pixel counts'):
        validate_pool_sum(cooler.Cooler(str(path)),[cooler.Cooler(s['matrix_path']+'::/resolutions/10000') for s in members],100000)


def test_depth_raw_only_conserves_draw_and_never_balances(prepared,monkeypatch):
    _,c=prepared
    c=dict(c,RUN_DEPTH_SENSITIVITY=True,RUN_HICREP=False,DOWNSAMPLING_REPEATS=1,KEEP_DOWNSAMPLED_MAPS=False)
    def forbidden(*a,**kw):raise AssertionError('raw-only thinning attempted balancing')
    monkeypatch.setattr(cooler,'balance_cooler',forbidden)
    before={s['sample_id']:file_hash(s['matrix_path']) for s in validate_metadata(c)['samples']}
    depth_sensitivity(c,'run')
    summary_path=Path(c['run_dir'])/'01_qc/depth_sensitivity/summary.tsv'
    first=pd.read_csv(summary_path,sep='\t')
    assert first.probability.between(0,1).all()
    assert (abs(first.realized_depth-first.target_expected_depth)/first.target_expected_depth<.02).all()
    assert not list((Path(c['run_dir'])/'01_qc/depth_sensitivity').rglob('*.cool'))
    depth_sensitivity(c,'run');pd.testing.assert_frame_equal(first,pd.read_csv(summary_path,sep='\t'))
    assert before=={s['sample_id']:file_hash(s['matrix_path']) for s in validate_metadata(c)['samples']}


def test_unchanged_copy_and_input_never_balanced(prepared,monkeypatch):
    root,c=prepared;c=dict(c,BALANCE_NEW_DERIVED_MAPS=True)
    source=cooler.Cooler(str(root/'C1.norm.mcool')+'::/resolutions/10000')
    copy=Path(c['run_dir'])/'02_pools/copy.cool';copy.parent.mkdir(parents=True,exist_ok=True)
    cooler.create_cooler(str(copy),source.bins()[:][['chrom','start','end']],source.pixels()[:],ordered=True)
    def forbidden(*a,**kw):raise AssertionError('Attempted balancing unchanged/imported counts')
    monkeypatch.setattr(cooler,'balance_cooler',forbidden)
    with pytest.raises(ValueError,match='unchanged'):balance_new_map(c,copy)
    with pytest.raises(ValueError,match='task-owned'):balance_new_map(c,root/'C1.norm.mcool')


def test_explicit_new_pools_balance_changed_owned_counts_only(prepared,monkeypatch):
    _,c=prepared;c=dict(c,RUN_POOL=True,POOL_MODE='create',BALANCE_NEW_DERIVED_MAPS=True,THREADS_BALANCE=1)
    samples=validate_metadata(c)['samples'];before={s['sample_id']:file_hash(s['matrix_path']) for s in samples}
    original=cooler.balance_cooler;seen=[]
    def guarded(clr,*args,**kwargs):
        path=Path(clr.filename).resolve()
        assert path.is_relative_to((Path(c['run_dir'])/'02_pools').resolve())
        assert clr.info['sum']>max(cooler.Cooler(s['matrix_path']+'::/resolutions/10000').info['sum'] for s in samples)
        seen.append(str(path));return original(clr,*args,**kwargs)
    monkeypatch.setattr(cooler,'balance_cooler',guarded)
    products=pools(c,'run')
    assert len(seen)==3
    registry=read_json(Path(c['run_dir'])/'02_pools/registry.json')['pools']
    assert all(set(row['uris'])=={'10000'} for row in registry)
    for row in registry:
        clr=cooler.Cooler(row['uris']['10000']);assert np.isfinite(clr.bins()['weight'][:]).any()
    assert before=={s['sample_id']:file_hash(s['matrix_path']) for s in samples}


def test_saddle_partition_ties_are_deterministic_without_pool(prepared):
    _,c=prepared;c=dict(c,COMPARTMENT_RESOLUTION=10000,SADDLE_PARTITION_SOURCE='replicate_rank_mean')
    samples=validate_metadata(c)['samples']
    bins=pd.read_csv(Path(c['run_dir'])/'00_metadata/masks/10000.bins.tsv',sep='\t')
    directory=Path(c['run_dir'])/'03_features/compartments';directory.mkdir(parents=True)
    for s in samples:
        table=bins.copy();table['oriented_score']=np.repeat(np.arange(60),10)
        table.to_csv(directory/(s['sample_id']+'.tsv'),sep='\t',index=False)
    first,_=saddle_partition(c,samples);second,_=saddle_partition(c,list(reversed(samples)))
    pd.testing.assert_frame_equal(first,second)
    assert set(first.partition)==set(range(50)) and first.partition.is_monotonic_increasing
    assert not (Path(c['run_dir'])/'02_pools').exists()


def test_cleanup_only_task_owned_and_requests_regeneration(prepared):
    root,c=prepared;c=dict(c,RUN_POOL=True,POOL_MODE='create',BALANCE_NEW_DERIVED_MAPS=True)
    publish_completion(c,'preflight','run',inspect_inputs(c,writing=True))
    publish_completion(c,'register','run',register(c,'run'))
    owned=Path(c['run_dir'])/'02_pools/new.cool';owned.parent.mkdir(parents=True,exist_ok=True);owned.write_bytes(b'task-owned')
    imported=root/'imported.cool';imported.write_bytes(b'imported')
    registry=Path(c['run_dir'])/'02_pools/registry.json'
    atomic_json(registry,{'pools':[{'pool_id':'new','task_owned':True,'uris':{'10000':str(owned)}},{'pool_id':'old','task_owned':False,'uris':{'10000':str(imported)}}]})
    publish_completion(c,'pool','run',[owned,registry])
    cleanup_derived(c)
    assert not owned.exists() and imported.read_bytes()==b'imported'
    assert complete(c,'pool','run')
    assert restore_required_pools(c,'loops') is True
    assert not record_path(c,'pool','run').exists()
