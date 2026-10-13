"""Required full-profile saddles and selective modules have distinct contracts."""
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from microc2compare import science
from microc2compare.contracts import ContractError,atomic_json,normalize,read_json
from microc2compare.analyses import resolution_set

ROOT=Path(__file__).resolve().parents[1]


def values(config):
    result={k:str(v).lower() if isinstance(v,bool) else str(v) for k,v in config.items()}
    result.update(ANALYSIS_PROFILE='practical_2_4kb',LOOP_CALLER='mustache',LOOP_RESOLUTION='2000',
                  LOOP_DISCOVERY_RESOLUTIONS='2000,4000',LOOP_QUANTIFICATION_RESOLUTIONS='2000,4000')
    return result


@pytest.mark.parametrize('profile',['practical_2_4kb','hsieh_inspired'])
@pytest.mark.parametrize('change',[{'RUN_COMPARTMENTS':'false'},{'RUN_AGGREGATES':'false'},{'AGGREGATE_MODULES':'apa,ada'},{'RUN_FINAL_REPORT':'false'}])
def test_full_profiles_reject_omitting_saddle_or_its_report(tmp_path,config,profile,change):
    env=values(config);env.update(ANALYSIS_PROFILE=profile,**change)
    with pytest.raises(ContractError,match='full fine-loop profile requires'):
        normalize(env,tmp_path)


def test_full_fine_profile_registers_compartment_grid_and_report_refresh_is_full(tmp_path,config):
    env=values(config);env.update(COMPARTMENT_RESOLUTION='100000')
    c=normalize(env,tmp_path)
    assert c['RUN_COMPARTMENTS'] and c['RUN_FINAL_REPORT'] and 'saddle' in c['AGGREGATE_MODULES'].split(',')
    assert 100000 in resolution_set(c)
    env.update(M2C_ANALYSES='report',AGGREGATE_MODULES='apa,ada')
    with pytest.raises(ContractError,match='full fine-loop profile requires'):normalize(env,tmp_path)


@pytest.mark.parametrize('selection',['qc,hicrep','saddle'])
def test_explicit_module_runs_keep_their_selected_prerequisites(tmp_path,config,selection):
    env=values(config);env.update(M2C_ANALYSES=selection,RUN_COMPARTMENTS='false',AGGREGATE_MODULES='apa,ada')
    c=normalize(env,tmp_path)
    if selection=='saddle':
        assert c['RUN_COMPARTMENTS'] and c['AGGREGATE_MODULES']=='saddle'
        assert not c['RUN_LOOPS'] and not c['RUN_INSULATION']
    else:
        assert not c['RUN_COMPARTMENTS'] and not c['RUN_AGGREGATES']


def test_missing_phasing_fails_preflight_before_contact_scan(tmp_path,config,monkeypatch):
    spec=importlib.util.spec_from_file_location('base_demo',ROOT/'scripts/generate_demo.py')
    generator=importlib.util.module_from_spec(spec);spec.loader.exec_module(generator)
    data=tmp_path/'data';generator.generate(data,development=True)
    registry=read_json(data/'references.json')
    registry.pop('phasing_track');registry.pop('phasing_track_sha256')
    atomic_json(data/'references.json',registry)
    env=values(config);env.update(M2C_ANALYSES='saddle',PROJECT='demo',REFERENCE_ID='demo',SAMPLESHEET='samples.csv',
                                 TECHNICAL_MEMBERSHIP='technical_membership.tsv',REFERENCE_REGISTRY='references.json',
                                 OUTDIR='results',SCRATCH_ROOT='scratch',TOTAL_CPUS='2',MAX_MEMORY_MB='16000')
    c=normalize(env,data)
    monkeypatch.setattr(science,'domain',lambda *args:pytest.fail('Contact scans must not precede required phasing validation'))
    with pytest.raises(ContractError,match='registered fixed phasing track'):science.inspect_inputs(c)


def test_no_eligible_saddle_partition_fails_instead_of_claiming_a_plot(tmp_path,config,monkeypatch):
    c=dict(config,run_dir=str(tmp_path),RUN_AGGREGATES=True,AGGREGATE_MODULES='saddle')
    monkeypatch.setattr(science,'validate_metadata',lambda c:{'samples':[]})
    monkeypatch.setattr(science,'saddle_partition',lambda *args:(pd.DataFrame({'partition':[-1,-1]}),np.array([0,0])))
    with pytest.raises(ContractError,match='Saddle has no eligible partition'):science._aggregates_single(c,'run')
