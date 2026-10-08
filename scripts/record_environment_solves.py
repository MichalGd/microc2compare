#!/usr/bin/env python3
"""Summarize actual dry-solve evidence without claiming Linux runtime tests."""
import argparse
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from microc2compare.contracts import atomic_json,digest,file_hash,read_json
p=argparse.ArgumentParser(description=__doc__)
for kind in ('core','hicrep','r','failed-r'):p.add_argument('--'+kind,required=True)
p.add_argument('--out',required=True);a=p.parse_args()
root=Path(__file__).resolve().parents[1];records={}
for kind,definition in [('core','environment.yml'),('hicrep','envs/hicrep.yml'),('r','envs/diff-r.yml')]:
    result=read_json(getattr(a,kind));packages=result.get('actions',{}).get('LINK',[])
    records[kind]={'conda_dry_solve_success':result.get('success',False),'runtime_verified_on_linux':False,'definition_sha256':file_hash(root/definition),
                   'package_count':len(packages),'solved_package_identity_digest':digest([{k:x.get(k) for k in ('name','version','build','url','sha256')} for x in packages]),
                   'principal_versions':{x['name']:x['version'] for x in packages if x['name'] in ('python','cooler','cooltools','bioframe','numpy','pandas','scipy','r-base','bioconductor-edger','bioconductor-limma','r-statmod','r-jsonlite')},
                   'pip_note':'HiCRep 0.2.6 installed and numerically exercised locally on Windows; Linux dry solve does not execute pip' if kind=='hicrep' else None}
failed=read_json(a.failed_r)
errors=[x['message'] for x in failed.get('log_history',[]) if x.get('level')=='error']
atomic_json(a.out,{'date':'2026-10-07','platform':'linux-64','solver':'micromamba 2.9.0 on Windows; dry cross-solve only','glibc_override':'2.17','records':records,
                  'rejected_r_candidate':{'r_base':'4.6','edgeR':'4.10.5','limma':'3.68.5','success':failed.get('success'), 'reason':errors},
                  'status':'dependency_solutions_available; full_Linux_runtime_and_exact_tested_locks_pending'})
