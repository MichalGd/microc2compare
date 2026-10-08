"""Resume refuses missing identities, wrong profiles and modified input files."""
import hashlib
import os
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def bash_path(path):
    if os.name=='nt': return '/'+path.drive[0].lower()+'/'+path.relative_to(path.anchor).as_posix()
    return str(path)


@pytest.mark.parametrize('problem,message', [
    ('missing_identity','original project.conf and input-sha256.txt'),
    ('wrong_profile','does not source this software profile'),
    ('changed_input','FAILED'),
])
def test_resume_rejects_without_regeneration_or_configuration_writes(tmp_path, problem, message):
    site=tmp_path/'software.conf'
    # Any accidental generator invocation fails because this executable is absent.
    site.write_text('CONDA_EXE="/must-not-be-called/conda"\nCORE_ENV_PREFIX="/unused"\n')
    destination=tmp_path/'existing'; destination.mkdir()
    matrix=destination/'identity-fixture.mcool'; matrix.write_bytes(b'Input identity fixture; not a scientific matrix')
    config=destination/'project.conf'
    config.write_text('source "'+bash_path(site)+'"\n')
    identity=destination/'input-sha256.txt'
    identity.write_text(hashlib.sha256(matrix.read_bytes()).hexdigest()+'  '+bash_path(matrix)+'\n')
    if problem=='missing_identity': identity.unlink()
    elif problem=='wrong_profile': config.write_text('source "/different/software.conf"\n')
    else: matrix.write_bytes(matrix.read_bytes()+b'modified')
    before={str(p):p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    bash=r'C:/Program Files/Git/bin/bash.exe' if os.name=='nt' else 'bash'
    result=subprocess.run([bash,str(ROOT/'scripts/validate_module_demo.sh'),
                           '--site',str(site),'--out',str(destination),'--resume'],capture_output=True,text=True)
    assert result.returncode!=0 and message in result.stdout+result.stderr
    assert 'must-not-be-called' not in result.stderr
    assert before=={str(p):p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
