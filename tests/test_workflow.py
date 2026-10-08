import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd
import pytest

from microc2compare.contracts import atomic_json,file_hash,normalize,read_json,validate_metadata
from microc2compare.science import context,hicrep,inspect_inputs,qc,register
from microc2compare.state import initialize

REPO=Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    spec=importlib.util.spec_from_file_location("demo_generator",REPO/"scripts/generate_demo.py")
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    root=tmp_path_factory.mktemp("demo")/"data"
    module.generate(root,development=True,raw_only=True)
    schema=read_json(REPO/"schemas/config.schema.json")
    values={k:str(v["default"]).lower() if isinstance(v["default"],bool) else str(v["default"]) for k,v in schema["properties"].items()}
    values.update(PROJECT="demo",ANALYSIS_SET_ID="synthetic_demo",RUN_ID="run001",REFERENCE_ID="demo",ENVIRONMENT_MODE="development",PYTHON_BOOTSTRAP=sys.executable,
                  REFERENCE_REGISTRY="references.json",INDEPENDENT_FEATURES_REGISTRY="",SAMPLESHEET="samples.csv",OUTDIR="results",SCRATCH_ROOT="scratch",
                  TOTAL_CPUS="2",MAX_MEMORY_MB="2000",THREADS_PREFLIGHT="1",MEMORY_MB_PREFLIGHT="1000",
                  RUN_COMPARTMENTS="false",RUN_INSULATION="false",RUN_LOOPS="false",RUN_AGGREGATES="false",RUN_DIFFERENTIAL="false")
    c=normalize(values,str(root))
    initialize(c,os.getpid(),"test")
    return root,c


def test_default_science_never_rebalances_or_copies(demo,monkeypatch):
    root,c=demo
    samples=validate_metadata(c)["samples"]
    before={s["sample_id"]:file_hash(s["matrix_path"]) for s in samples}
    def forbidden(*args,**kwargs):raise AssertionError("Default run attempted new matrix construction/balancing")
    for name in ("balance_cooler","merge_coolers","zoomify_cooler","coarsen_cooler","create_cooler"):
        monkeypatch.setattr(cooler,name,forbidden)
    inspect_inputs(c,writing=True);register(c,"run")
    for sample in samples:qc(c,sample["sample_id"])
    hicrep(c,"run")
    assert before=={s["sample_id"]:file_hash(s["matrix_path"]) for s in samples}
    assert not list(Path(c["run_dir"]).rglob("*.cool")) and not list(Path(c["run_dir"]).rglob("*.mcool"))
    summary=pd.read_csv(Path(c["run_dir"])/"01_qc/hicrep/summary.tsv",sep="\t")
    assert len(summary)==15 and summary.weighted_scc.between(-1,1).all()


def test_missing_weights_fail_without_repair(demo,monkeypatch):
    root,c=demo
    changed=dict(c,WEIGHT_COLUMN="missing_weight")
    def forbidden(*args,**kwargs):raise AssertionError("Unexpected balancing repair")
    monkeypatch.setattr(cooler,"balance_cooler",forbidden)
    with pytest.raises(ValueError,match="missing_weight"):inspect_inputs(changed)


def test_optional_switch_alone_does_not_create_tasks(demo):
    _,c=demo
    from microc2compare.contracts import enabled
    c=dict(c,BALANCE_NEW_DERIVED_MAPS=True)
    assert not enabled(c,"pool") and not enabled(c,"depth_sensitivity")


def test_background_child_failure_is_propagated():
    bash=shutil.which("bash")
    if os.name=="nt":bash=r"C:\Program Files\Git\bin\bash.exe"
    if not bash or not Path(bash).exists():pytest.skip("Bash unavailable")
    text='source scripts/lib/common.sh; source scripts/lib/parallel.sh; worker(){ [[ "$1" != fail ]]; }; if batch_run 2 worker ok fail; then exit 0; else exit 9; fi'
    result=subprocess.run([bash,"-c",text],cwd=REPO,capture_output=True,text=True)
    assert result.returncode==9


def test_raw_only_hicrep_can_read_maps_without_weights(tmp_path,monkeypatch):
    spec=importlib.util.spec_from_file_location('generator',REPO/'scripts/generate_demo.py')
    generator=importlib.util.module_from_spec(spec);spec.loader.exec_module(generator)
    root=tmp_path/'data';generator.generate(root,development=True,raw_only=True)
    for path in root.glob('*.mcool'):
        with h5py.File(path,'r+') as handle:
            for group in handle['resolutions'].values():del group['bins/weight']
    schema=read_json(REPO/'schemas/config.schema.json')
    values={k:str(v['default']).lower() if isinstance(v['default'],bool) else str(v['default']) for k,v in schema['properties'].items()}
    values.update(PROJECT='demo',ANALYSIS_SET_ID='synthetic_demo',RUN_ID='raw',REFERENCE_ID='demo',ENVIRONMENT_MODE='development',PYTHON_BOOTSTRAP=sys.executable,
                  REFERENCE_REGISTRY='references.json',SAMPLESHEET='samples.csv',OUTDIR='results',SCRATCH_ROOT='scratch',TOTAL_CPUS='2',MAX_MEMORY_MB='2000',
                  RUN_QC='false',RUN_COMPARTMENTS='false',RUN_INSULATION='false',RUN_LOOPS='false',RUN_AGGREGATES='false',RUN_DIFFERENTIAL='false')
    c=normalize(values,str(root));initialize(c,os.getpid(),'test')
    before={p.name:file_hash(p) for p in root.glob('*.mcool')}
    monkeypatch.setattr(cooler,'balance_cooler',lambda *a,**kw:pytest.fail('raw-only HiCRep attempted balancing'))
    inspect_inputs(c,writing=True);register(c,'run');hicrep(c,'run')
    assert before=={p.name:file_hash(p) for p in root.glob('*.mcool')}
