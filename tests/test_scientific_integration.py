"""Actual core cooltools APIs, with construction/balancing forbidden afterward."""
import importlib.util
import os
import sys
from pathlib import Path

import cooler
import numpy as np
import pandas as pd
import pytest

from microc2compare.contracts import atomic_json,file_hash,normalize,read_json,validate_metadata
from microc2compare.science import IMPLEMENTATIONS,saddle_partition
from microc2compare.state import initialize,publish_completion,scopes

REPO=Path(__file__).resolve().parents[1]


@pytest.mark.scientific
@pytest.mark.parametrize('mixed_scope',[False,True])
def test_actual_cooltools_masked_read_only_adapters(tmp_path,monkeypatch,mixed_scope):
    pytest.importorskip('cooltools',reason='Full Linux cooltools runtime remains a separate gate')
    # These stages run in core; HiCRep is intentionally in a separate prefix.
    spec=importlib.util.spec_from_file_location('generator',REPO/'scripts'/('generate_common_chromosome_demo.py' if mixed_scope else 'generate_demo.py'))
    generator=importlib.util.module_from_spec(spec);spec.loader.exec_module(generator)
    root=tmp_path/'data';generator.generate(root,development=True)
    blacklist=root/'blacklist.bed';blacklist.write_text('chrX\t100000\t110000\n' if mixed_scope else 'chrDemo\t2500000\t2510000\n')
    reference=read_json(root/'references.json');reference.update(blacklist=str(blacklist),blacklist_sha256=file_hash(blacklist));atomic_json(root/'references.json',reference)
    schema=read_json(REPO/'schemas/config.schema.json')
    values={k:str(v['default']).lower() if isinstance(v['default'],bool) else str(v['default']) for k,v in schema['properties'].items()}
    values.update(PROJECT='demo',ANALYSIS_SET_ID='synthetic_demo',RUN_ID='run001',REFERENCE_ID='demo',ENVIRONMENT_MODE='development',PYTHON_BOOTSTRAP=sys.executable,
                  REFERENCE_REGISTRY='references.json',SAMPLESHEET='samples.csv',OUTDIR='results',SCRATCH_ROOT='scratch',TOTAL_CPUS='2',MAX_MEMORY_MB='2000')
    for key in values:
        if key.startswith('THREADS_'):values[key]='1'
        elif key.startswith('MEMORY_MB_'):values[key]='1000'
    c=normalize(values,str(root));initialize(c,os.getpid(),'test')
    samples=validate_metadata(c)['samples'];before={s['sample_id']:file_hash(s['matrix_path']) for s in samples}
    def forbidden(*a,**kw):raise AssertionError('Core actual-tool run tried matrix construction or balancing')
    for name in ('balance_cooler','merge_coolers','zoomify_cooler','coarsen_cooler','create_cooler'):monkeypatch.setattr(cooler,name,forbidden)
    for stage in ('preflight','register','qc','compartments','insulation','loops','consensus','tads','quantify','aggregates','differential'):
        for scope in scopes(c,stage):
            products=IMPLEMENTATIONS[stage](c,scope)
            publish_completion(c,stage,scope,products)
    partition,_=saddle_partition(c,samples);again,_=saddle_partition(c,list(reversed(samples)))
    pd.testing.assert_frame_equal(partition,again)
    assert set(partition.partition)==set(range(50))|{-1}
    for s in samples:
        table=pd.read_csv(Path(c['run_dir'])/f"03_features/compartments/{s['sample_id']}.tsv",sep='\t')
        phase=pd.read_csv(root/'phasing.bedgraph',sep='\t',header=None).iloc[:,3]
        good=np.isfinite(table.oriented_score)
        assert np.corrcoef(table.loc[good,'oriented_score'],phase[good])[0,1]>.1
        assert pd.isna(table.loc[21 if mixed_scope else 25,'oriented_score'])
        insulation=pd.read_csv(Path(c['run_dir'])/f"03_features/insulation/{s['sample_id']}.tsv",sep='\t')
        assert insulation.log2_insulation_score_100000.notna().sum()>400
    assert before=={s['sample_id']:file_hash(s['matrix_path']) for s in samples}
    assert not list(Path(c['run_dir']).rglob('*.cool')) and not list(Path(c['run_dir']).rglob('*.mcool'))
