"""Real state/validator CLIs under a stdlib-only bootstrap (-I -S).

Completion fixtures test coordination and hashes, not contact-map numerics.
The scientific environment must never be needed to inspect their identities.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from microc2compare.contracts import (SAMPLE_FIELDS, atomic_json, normalize,
                                    write_table)
from microc2compare.state import initialize, publish_completion

ROOT = Path(__file__).resolve().parents[1]
BOOTSTRAP_WRAPPER = '''
import importlib.abc
import runpy
import sys
scientific = {'cooler','numpy','pandas','h5py','cooltools','bioframe'}
assert not scientific.intersection(sys.modules)
class StandardLibraryBootstrap(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in scientific:
            raise ModuleNotFoundError('Bootstrap has no scientific package: '+fullname)
sys.meta_path.insert(0, StandardLibraryBootstrap())
sys.argv = sys.argv[1:]
runpy.run_path(sys.argv[0], run_name='__main__')
'''


@pytest.fixture
def registered_state(tmp_path, config):
    matrix = tmp_path/'input-identity-fixture.txt'
    matrix.write_text('State-contract input fixture; not a scientific contact matrix\n')
    provenance = tmp_path/'provenance.json'
    atomic_json(provenance, {})
    sample = dict(sample_id='C1', project='demo', condition='control', assay='microc',
                  biological_replicate='1', subject_id='donor1', reference_id='demo',
                  matrix_path=str(matrix), provenance_path=str(provenance))
    write_table(tmp_path/'samples.csv', [sample], SAMPLE_FIELDS)
    member = dict(sample_id='C1', library_id='lib1', source_id='source1',
                  source_kind='independent_library', dedup_group_id='lib1',
                  dedup_scope='joint_library', source_checksum='a'*64)
    write_table(tmp_path/'technical_membership.tsv', [member], list(member))
    atomic_json(tmp_path/'references.json', {'reference_id':'demo'})
    values = {k: str(v).lower() if isinstance(v, bool) else str(v) for k,v in config.items()}
    for key in values:
        if key.startswith('RUN_'): values[key]='false'
        elif key.startswith('THREADS_'): values[key]='1'
        elif key.startswith('MEMORY_MB_'): values[key]='1000'
    values.update(PROJECT='demo', ANALYSIS_SET_ID='bootstrap', RUN_ID='run001',
                  REFERENCE_ID='demo', REFERENCE_REGISTRY='references.json',
                  SAMPLESHEET='samples.csv', OUTDIR='results', SCRATCH_ROOT='scratch',
                  TOTAL_CPUS='2', MAX_MEMORY_MB='2000', ENVIRONMENT_MODE='development',
                  PYTHON_BOOTSTRAP=sys.executable)
    c = normalize(values, tmp_path)
    initialize(c, os.getpid(), 'bootstrap-test')
    marker = Path(c['run_dir'])/'00_metadata/coordination-fixture.json'
    atomic_json(marker, {'purpose':'State dependency and output-hash regression only'})
    publish_completion(c, 'preflight', 'run', [marker])
    publish_completion(c, 'register', 'run', [marker])
    return c, Path(c['run_dir'])/'00_metadata/config.json'


@pytest.mark.parametrize('command,arguments', [
    ('state.py', ['prerequisites','--stage','qc']),
    ('state.py', ['summarize']),
    ('state.py', ['failed','--message','intentional test status']),
    ('validate_outputs.py', ['--until-stage','register']),
])
def test_completion_cli_without_scientific_site_packages(registered_state, command, arguments):
    c, snapshot = registered_state
    option = '--config-json' if command=='validate_outputs.py' else '--config'
    result = subprocess.run([sys.executable, '-I', '-S', '-c', BOOTSTRAP_WRAPPER, str(ROOT/'scripts'/command),
                             *arguments, option, str(snapshot)], capture_output=True, text=True)
    assert result.returncode==0, result.stdout+result.stderr
    if command=='validate_outputs.py':
        assert 'register/run: valid' in result.stdout
    elif arguments[0] in ('summarize','failed'):
        summary=json.loads((Path(c['run_dir'])/'run_summary.json').read_text())
        assert len(summary['tasks'])==2 and all(r['currently_valid'] for r in summary['tasks'])
        assert summary['status']==('failed' if arguments[0]=='failed' else 'completed_requested_stages')
