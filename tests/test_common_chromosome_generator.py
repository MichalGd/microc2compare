"""Validate the server acceptance fixture with real maps and production preflight."""
import importlib.util
import os
import sys
from pathlib import Path

import cooler
import pandas as pd
import pytest

from microc2compare.contracts import file_hash, normalize, read_json, validate_metadata
from microc2compare.science import inspect_inputs, qc, register
from microc2compare.state import initialize


def test_new_fixture_preflight_common_xy_qc_and_no_overwrite(tmp_path):
    repo=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('common_generator',repo/'scripts/generate_common_chromosome_demo.py')
    generator=importlib.util.module_from_spec(spec);spec.loader.exec_module(generator)
    root=tmp_path/'data'
    generator.generate(root,development=True)
    schema=read_json(repo/'schemas/config.schema.json')
    values={k:str(v['default']).lower() if isinstance(v['default'],bool) else str(v['default']) for k,v in schema['properties'].items()}
    values.update(PROJECT='demo',ANALYSIS_SET_ID='common_chromosomes',RUN_ID='check',REFERENCE_ID='demo',
                  ENVIRONMENT_MODE='development',PYTHON_BOOTSTRAP=sys.executable,REFERENCE_REGISTRY='references.json',
                  SAMPLESHEET='samples.csv',OUTDIR='results',SCRATCH_ROOT='scratch',TOTAL_CPUS='2',MAX_MEMORY_MB='2000',
                  RUN_COMPARTMENTS='false',RUN_INSULATION='false',RUN_LOOPS='false',RUN_AGGREGATES='false',RUN_DIFFERENTIAL='false')
    values.update(CHROMOSOME_SCOPE='autosomes_plus_sex',SEX_CHROMOSOMES='chrX,chrY')
    c=normalize(values,str(root));initialize(c,os.getpid(),'test')
    samples=validate_metadata(c)['samples'];before={s['sample_id']:file_hash(s['matrix_path']) for s in samples}
    inspect_inputs(c,writing=True);register(c,'run')
    for sample in samples:
        qc(c,sample['sample_id'])
        raw=cooler.Cooler(sample['matrix_path']+'::/resolutions/10000')
        assert raw.chromnames==(['scaffold','chrY','chr1','chrX'] if sample['sample_id'].startswith('C') else ['chr1','chrX','chrY'])
        ps=pd.read_csv(Path(c['run_dir'])/f"01_qc/ps/{sample['sample_id']}.10000.tsv",sep='\t')
        assert set(ps.region1)=={'chr1','chrX','chrY'}
    with pytest.raises(ValueError,match='empty'):generator.generate(root,development=True)
    assert before=={s['sample_id']:file_hash(s['matrix_path']) for s in samples}
