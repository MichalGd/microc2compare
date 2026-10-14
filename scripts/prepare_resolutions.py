#!/usr/bin/env python3
"""Explicit preparation of missing integer coarsenings; imported maps stay read-only."""
import argparse
import copy
import json
import os
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import cooler
import numpy as np
from microc2compare.contracts import IDENTIFIER,atomic_json,file_hash,read_json,read_table,require,resolve,write_table,load_structured
from microc2compare.matrix import matrix_uri
from microc2compare.resolutions import sample_matrix


def prepare(samplesheet,out,resolutions,threads=4,max_iters=1000,matrix_registry=''):
    os.umask(0o077)
    root=Path(out).resolve();require(not root.exists(),'Preparation destination exists; choose a new directory')
    samples=read_table(samplesheet,('sample_id','matrix_path','provenance_path','reference_id'))
    require(samples and len({s['sample_id'] for s in samples})==len(samples),'Invalid/duplicate biological samples')
    require(all(IDENTIFIER.fullmatch(s['sample_id']) for s in samples),'Unsafe sample identifier')
    require(resolutions and len(resolutions)==len(set(resolutions)) and all(r>0 for r in resolutions),'Invalid/duplicate requested resolutions')
    require(len({s['reference_id'] for s in samples})==1,'Mixed reference identities')
    if matrix_registry:
        require(load_structured(matrix_registry).get('reference_id')==samples[0]['reference_id'],'Matrix registry reference mismatch')
    root.mkdir(parents=True);registry={'reference_id':samples[0]['reference_id'],'samples':{}}
    for sample in samples:
        source=resolve(sample['matrix_path'],Path(samplesheet).parent)
        original_provenance=resolve(sample['provenance_path'],Path(samplesheet).parent)
        prov=copy.deepcopy(read_json(original_provenance));before=file_hash(source.split('::')[0])
        configured={'REFERENCE_ID':sample['reference_id'],'MATRIX_REGISTRY':matrix_registry}
        registered=load_structured(matrix_registry) if matrix_registry else {}
        other_sources=[resolve(value,Path(matrix_registry).parent) for value in registered.get('samples',{}).get(sample['sample_id'],{}).values()] if matrix_registry else []
        identities={source.split('::')[0]:before}
        for path in other_sources:
            file=path.split('::')[0]
            if file not in identities:identities[file]=file_hash(file)
        registry['samples'][sample['sample_id']]={}
        for res in resolutions:
            try:
                uri=matrix_uri(sample_matrix(configured,sample | {'matrix_path':source},res),res)
                registry['samples'][sample['sample_id']][str(res)]=uri
                print(f"{sample['sample_id']} {res}: reuse {uri}",flush=True);continue
            except ValueError:pass
            available=[path+'::'+g for path in identities for g in cooler.fileops.list_coolers(path)]
            eligible=[u for u in available if cooler.Cooler(u).binsize and cooler.Cooler(u).binsize<res and res%cooler.Cooler(u).binsize==0]
            require(eligible,f"{sample['sample_id']}: no exact finer source for {res}; rebuild upstream from suitable contact pairs")
            base=min(eligible,key=lambda u:cooler.Cooler(u).binsize)
            original=cooler.Cooler(base);factor=int(res//original.binsize)
            require(original.pixels()['count'][0:0].dtype.kind in 'iu' and original.info.get('storage-mode')=='symmetric-upper','Preparation requires symmetric-upper raw integer counts')
            target=root/f"{sample['sample_id']}.{res}.cool"
            print(f"{sample['sample_id']} {res}: coarsen {base}, factor={factor}",flush=True)
            cooler.coarsen_cooler(base,str(target),factor=factor,chunksize=1000000,nproc=threads)
            derived=cooler.Cooler(str(target))
            require(int(derived.info['sum'])==int(original.info['sum']),'Raw contact conservation failed')
            with original.open('r') as h:
                require('weight' in h['bins'],'Source balancing recipe unavailable')
                attrs=dict(h['bins/weight'].attrs)
            keys=('cis_only','ignore_diags','mad_max','min_count','min_nnz','tol')
            require(all(k in attrs for k in keys),'Source stored-weight recipe incomplete; supply corrected upstream metadata')
            require(not bool(attrs.get('divisive_weights',False)) and not bool(attrs.get('trans_only',False)),'Unsupported source balancing convention')
            recipe={k:attrs[k].item() if isinstance(attrs[k],np.generic) else attrs[k] for k in keys}
            from multiprocessing import Pool
            with Pool(threads) as workers:
                _,stats=cooler.balance_cooler(derived,store=True,map=workers.map,chunksize=1000000,max_iters=max_iters,**recipe)
            require(bool(np.asarray(stats['converged']).all()),f'{target}: new balancing failed; prepared inputs are not ready')
            prov.setdefault('balancing',{})[str(res)]={'weight_column':'weight','convention':'multiplicative','converged':True,'recipe':recipe}
            prov.setdefault('derived_resolutions',{})[str(res)]={'source_uri':base,'source_file_sha256':identities[base.split('::')[0]],'factor':factor,'native_contact_count':int(derived.info['sum']),'max_iters':max_iters,'generated_by':f"cooler {cooler.__version__}"}
            registry['samples'][sample['sample_id']][str(res)]=str(target)
        require(all(file_hash(path)==sha for path,sha in identities.items()),'A reused matrix changed during preparation')
        provenance=root/(sample['sample_id']+'.provenance.json');atomic_json(provenance,prov)
        sample['matrix_path']=source;sample['provenance_path']=str(provenance)
    atomic_json(root/'matrix_registry.json',registry)
    write_table(root/'samples.tsv',samples,list(samples[0]))
    atomic_json(root/'preparation.json',{'status':'prepared; run downstream preflight and fine-resolution pilot','resolutions':resolutions,'samples':len(samples),'inputs_modified':False})
    print(f"Set SAMPLESHEET={root/'samples.tsv'} and MATRIX_REGISTRY={root/'matrix_registry.json'}",flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--samples',required=True);p.add_argument('--out',required=True);p.add_argument('--resolutions',default='2000,4000');p.add_argument('--threads',type=int,default=4);p.add_argument('--max-iters',type=int,default=1000)
    a=p.parse_args();require(a.threads>0 and a.max_iters>0,'Threads and iterations must be positive');prepare(a.samples,a.out,[int(x) for x in a.resolutions.split(',')],a.threads,a.max_iters)
