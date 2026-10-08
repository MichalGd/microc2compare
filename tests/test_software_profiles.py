import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
BASH=Path('C:/Program Files/Git/bin/bash.exe') if os.name=='nt' else Path(shutil.which('bash') or '/missing/bash')


def shell_path(path):
    path=Path(path).resolve()
    if os.name=='nt':
        return subprocess.check_output([str(BASH),'-c','cygpath -u "$1"','bash',str(path)],text=True).strip()
    return str(path)


def module(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/(name+'.py'))
    value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value);return value


@pytest.fixture
def installation(tmp_path):
    if not BASH.is_file():pytest.skip('Bash unavailable')
    base=tmp_path/'shared software'
    archive=tmp_path/'release.zip'
    subprocess.run([sys.executable,str(ROOT/'scripts/package_release.py'),'--out',str(archive)],check=True,capture_output=True)
    source=module('install_release').install(archive,hashlib.sha256(archive.read_bytes()).hexdigest(),base,'test')
    site=base/'site';site.mkdir();site.chmod(0o755)
    profile=site/'software profile.conf'
    profile.write_text('ENVIRONMENT_MODE="development"\nPYTHON_BOOTSTRAP='+json.dumps(Path(sys.executable).as_posix())+'\nCORE_ENV_PREFIX='+json.dumps((base/'checked core').as_posix())+'\nDIFF_R_ENV_PREFIX='+json.dumps((base/'checked r').as_posix())+'\n',newline='\n')
    profile.chmod(0o644)
    publisher=module('publish_shared_command')
    publisher.publish(source,base,site_profile=profile)
    launcher=base/'bin'/('microc2compare-'+json.loads((source/'RELEASE_MANIFEST.json').read_text())['version']+'-test')
    env={k:v for k,v in os.environ.items() if not k.startswith('M2C_')}
    env['PYTHONDONTWRITEBYTECODE']='1'
    return base,source,profile,launcher,publisher,env


def load(installation,project,override=None):
    base,source,profile,_,_,env=installation
    env=dict(env,M2C_DEFAULT_SOFTWARE_PROFILE=profile.as_posix(),
             M2C_DEFAULT_SOFTWARE_PROFILE_SHA256=hashlib.sha256(profile.read_bytes()).hexdigest())
    if override is not None:env['M2C_SOFTWARE_PROFILE']=override.as_posix()
    code='source "$1/scripts/lib/common.sh"; load_config "$2"; cat "$(shell_path "$CONFIG_JSON")"; rm -- "$(shell_path "$CONFIG_JSON")"'
    return subprocess.run([str(BASH),'-c',code,'bash',shell_path(source),shell_path(project)],env=env,capture_output=True,text=True)


def test_defaults_project_analysis_and_software_precedence(installation,tmp_path):
    project=tmp_path/'private project'/'project.conf';project.parent.mkdir()
    project.write_text('PROJECT="my_project"\nRUN_ID="my_run"\nTOTAL_CPUS="2"\nMAX_MEMORY_MB="2000"\nRUN_TADS="false"\nOUTDIR="my-results"\nCORE_ENV_PREFIX="/placeholder"\nPYTHON_BOOTSTRAP="/missing/python"\n')
    result=load(installation,project)
    assert result.returncode==0,result.stderr
    config=json.loads(result.stdout)
    assert config['PROJECT']=='my_project' and config['RUN_ID']=='my_run'
    assert config['TOTAL_CPUS']==2 and config['RUN_TADS'] is False
    assert config['OUTDIR']==str(project.parent/'my-results')
    assert config['CORE_ENV_PREFIX']==str(installation[0]/'checked core') and config['DIFF_R_ENV_PREFIX']==str(installation[0]/'checked r')
    assert config['_software_profile']['selection']=='launcher_default'
    assert config['_software_profile']['sha256']==hashlib.sha256(installation[2].read_bytes()).hexdigest()


def test_profile_cannot_replace_project_analysis_values(installation,tmp_path):
    override=tmp_path/'override.conf'
    override.write_text(installation[2].read_text()+'printf "PROJECT=stdout_injection\\n"\nPROJECT="unwanted"\nTOTAL_CPUS=99\nOUTDIR="/unwanted"\nRUN_QC=false\nCORE_ENV_PREFIX='+json.dumps((tmp_path/'override core').as_posix())+'\n')
    project=tmp_path/'project.conf';project.write_text('PROJECT="mine"\nTOTAL_CPUS=2\nOUTDIR="results"\n')
    result=load(installation,project,override)
    assert result.returncode==0,result.stderr
    config=json.loads(result.stdout)
    assert config['PROJECT']=='mine' and config['TOTAL_CPUS']==2 and config['RUN_QC'] is True
    assert config['OUTDIR']==str(tmp_path/'results') and config['CORE_ENV_PREFIX']==str(tmp_path/'override core')
    assert config['_software_profile']['selection']=='explicit_override'


def test_actual_launcher_list_uses_bound_bootstrap_and_cli_override(installation,tmp_path):
    _,_,profile,launcher,_,env=installation
    # Hide python3 from PATH: the bound absolute bootstrap must be used.
    env['PYTHON_BOOTSTRAP']='/missing/python'
    result=subprocess.run([str(BASH),shell_path(launcher),'--list-analyses'],env=env,capture_output=True,text=True)
    assert result.returncode==0 and 'tads' in result.stdout,result.stderr
    override=tmp_path/'override with spaces.conf'
    override.write_bytes(profile.read_bytes())
    # Options after --list-analyses are parsed too; relative paths use invocation cwd.
    result=subprocess.run([str(BASH),shell_path(launcher),'--list-analyses','--software-profile',override.name],cwd=tmp_path,env=env,capture_output=True,text=True)
    assert result.returncode==0 and 'compartments' in result.stdout,result.stderr


def test_changed_default_profile_fails_and_explicit_override_remains_available(installation,tmp_path):
    _,_,profile,launcher,_,env=installation
    original=profile.read_bytes();profile.write_bytes(original+b'# changed\n')
    result=subprocess.run([str(BASH),shell_path(launcher),'--list-analyses'],env=env,capture_output=True,text=True)
    assert result.returncode!=0 and 'Pinned software profile changed' in result.stderr
    override=tmp_path/'reviewed.conf';override.write_bytes(original)
    result=subprocess.run([str(BASH),shell_path(launcher),'--software-profile',str(override),'--list-analyses'],env=env,capture_output=True,text=True)
    assert result.returncode==0,result.stderr


@pytest.mark.parametrize('content,message',[(None,'not readable'),('false\n','could not be loaded')])
def test_missing_or_failed_override_stops_configuration(installation,tmp_path,content,message):
    profile=tmp_path/'override.conf'
    if content is not None:profile.write_text(content)
    project=tmp_path/'project.conf';project.write_text('PROJECT="mine"\n')
    result=load(installation,project,profile)
    assert result.returncode!=0 and message in result.stderr
    assert not (tmp_path/'results').exists()


def test_profile_values_remain_literal(installation,tmp_path):
    profile=tmp_path/'literal.conf';sentinel=tmp_path/'should-not-exist'
    literal='$(touch '+shell_path(sentinel)+')'
    profile.write_text(installation[2].read_text()+"CORE_ENV_PREFIX='"+literal+"'\n")
    project=tmp_path/'project.conf';project.write_text('PROJECT="mine"\n')
    result=load(installation,project,profile)
    assert result.returncode==0,result.stderr
    assert not sentinel.exists()
    assert '$(touch' in json.loads(result.stdout)['CORE_ENV_PREFIX']


def test_publisher_refuses_external_profile_and_rebinding(installation,tmp_path):
    base,source,profile,launcher,publisher,_=installation
    outside=tmp_path/'outside.conf';outside.write_bytes(profile.read_bytes())
    with pytest.raises(ValueError,match='below the shared site'):publisher.publish(source,base,site_profile=outside)
    second=profile.parent/'second.conf';second.write_bytes(profile.read_bytes());second.chmod(0o644)
    before=launcher.read_bytes()
    with pytest.raises(ValueError,match='refusing overwrite'):publisher.publish(source,base,site_profile=second)
    assert launcher.read_bytes()==before


def test_unbound_direct_invocation_keeps_project_software_and_clears_stale_record(installation,tmp_path):
    _,source,_,_,_,env=installation
    project=tmp_path/'project.conf'
    project.write_text('ENVIRONMENT_MODE="development"\nPYTHON_BOOTSTRAP='+json.dumps(Path(sys.executable).as_posix())+'\n')
    env=dict(env,M2C_LOADED_SOFTWARE_PROFILE='/stale/profile',
             M2C_LOADED_SOFTWARE_PROFILE_SHA256='stale',M2C_SOFTWARE_PROFILE_ORIGIN='launcher_default')
    code='source "$1/scripts/lib/common.sh"; load_config "$2"; cat "$(shell_path "$CONFIG_JSON")"; rm -- "$(shell_path "$CONFIG_JSON")"'
    result=subprocess.run([str(BASH),'-c',code,'bash',shell_path(source),shell_path(project)],env=env,capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    config=json.loads(result.stdout)
    assert config['ENVIRONMENT_MODE']=='development' and '_software_profile' not in config


def test_plain_config_actual_qc_run_uses_launcher_profile(installation,tmp_path):
    base,source,profile,launcher,_,env=installation
    data=tmp_path/'private demo'
    module('generate_demo').generate(data,development=True,raw_only=True)
    original=data/'project.conf'
    text=original.read_text()
    # Remove all software settings: the actual launcher must supply them.
    software={'PYTHON_BOOTSTRAP','CONDA_EXE','ENVIRONMENT_MODE','CORE_ENV_PREFIX','HICREP_ENV_PREFIX','DIFF_R_ENV_PREFIX','CORE_LOCK','HICREP_LOCK','DIFF_R_LOCK','R_SCRIPT'}
    text='\n'.join(line for line in text.splitlines() if line.split('=',1)[0] not in software)
    text+='\nTOTAL_CPUS="2"\nMAX_MEMORY_MB="2000"\nMEMORY_MB_PREFLIGHT="1000"\nMEMORY_MB_QC="1000"\nMEMORY_MB_REPORT="1000"\n'
    original.write_text(text+'\n',newline='\n')
    before={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in data.glob('*.mcool')}
    result=subprocess.run([str(BASH),shell_path(launcher),'-c',shell_path(original),'--analysis','qc'],env=env,capture_output=True,text=True,timeout=600)
    assert result.returncode==0,result.stderr
    report=data/'results/demo/synthetic_demo/run001/08_reports/index.html'
    assert report.is_file()
    config=json.loads((report.parents[1]/'00_metadata/config.json').read_text())
    assert config['_software_profile']['selection']=='launcher_default'
    assert Path(config['PYTHON_BOOTSTRAP'])==Path(sys.executable)
    assert before=={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in data.glob('*.mcool')}
    assert not list((data/'results').rglob('*.mcool'))


@pytest.mark.skipif(os.name=='nt',reason='Unix write permission bits are not Windows ACLs')
def test_publisher_rejects_group_writable_profile(installation):
    base,source,profile,_,publisher,_=installation
    profile.chmod(0o664)
    with pytest.raises(ValueError,match='writable by non-maintainers'):
        publisher.publish(source,base,site_profile=profile)
