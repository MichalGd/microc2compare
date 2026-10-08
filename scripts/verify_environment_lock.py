#!/usr/bin/env python3
"""Dependency-free verification of registered immutable prefix inventories."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from microc2compare.contracts import digest,read_json,require
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True);p.add_argument('--kind',choices=('core','hicrep','diff_r'),required=True);a=p.parse_args()
c=read_json(a.config)
if c['ENVIRONMENT_MODE']=='locked':
    prefix_key,lock_key={'core':('CORE_ENV_PREFIX','CORE_LOCK'),'hicrep':('HICREP_ENV_PREFIX','HICREP_LOCK'),'diff_r':('DIFF_R_ENV_PREFIX','DIFF_R_LOCK')}[a.kind]
    lock=c[lock_key]
    require(lock and Path(lock+'.inventory.json').is_file(),f'Locked {a.kind} requires an exact lock inventory')
    actual=subprocess.run([c['CONDA_EXE'],'list','-p',c[prefix_key],'--json'],capture_output=True,text=True,check=True)
    require(digest(json.loads(actual.stdout))==read_json(lock+'.inventory.json')['package_digest'],f'{a.kind} environment drift from registered exact lock')
