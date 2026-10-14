#!/usr/bin/env python3
"""Check pinned numeric APIs; does not install, rebalance or run user analyses."""
import argparse
import importlib.metadata
import json
import numpy as np
p=argparse.ArgumentParser()
p.add_argument('--mustache',action=argparse.BooleanOptionalAction,default=False)
p.add_argument('--chromosight',action=argparse.BooleanOptionalAction,default=False)
a=p.parse_args();versions={}
if a.mustache:
    from mustache import mustache
    assert all(callable(getattr(mustache,n)) for n in ('normalize_sparse','mustache'))
    assert importlib.metadata.version('mustache-hic')=='1.3.3'
    versions['mustache-hic']='1.3.3'
if a.chromosight:
    from chromosight.utils.detection import normxcorr2
    from importlib.resources import files
    assert importlib.metadata.version('chromosight')=='1.6.3'
    kernel=np.loadtxt(str(files('chromosight.kernels').joinpath('artificial_template_loops_type1.txt')))
    signal=np.zeros((kernel.shape[0]+4,kernel.shape[1]+4));signal[2:-2,2:-2]=kernel
    result,_=normxcorr2(signal,kernel,pval=False)
    assert np.nanmax(result)>.99
    versions['chromosight']='1.6.3'
print(json.dumps({'versions':versions,'status':'numeric_API_probe_only; Linux and real-data acceptance remain separate'}))
