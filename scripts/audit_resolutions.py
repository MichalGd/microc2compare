#!/usr/bin/env python3
"""Read-only resolution inventory; missing grids are reported, never generated."""
import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import cooler
import numpy as np
import pandas as pd
from microc2compare.contracts import atomic_json,read_table,resolve,require,read_json,load_structured,file_hash
from microc2compare.matrix import matrix_uri,reference,native_indices,analysis_bins
from microc2compare.resolutions import sample_matrix


def audit(samplesheet,reference_registry,out,resolutions,matrix_registry='',scope='autosomes',sex='',chromosomes=''):
    destination=Path(out).resolve();require(not destination.exists(),'Choose a NEW audit directory')
    registry=load_structured(reference_registry)
    c={k:v['default'] for k,v in read_json(Path(__file__).resolve().parents[1]/'schemas/config.schema.json')['properties'].items()}
    c.update(REFERENCE_ID=registry['reference_id'],REFERENCE_REGISTRY=str(Path(reference_registry).resolve()),RUN_COMPARTMENTS=False,
             MATRIX_REGISTRY=str(Path(matrix_registry).resolve()) if matrix_registry else '',CHROMOSOME_SCOPE=scope,SEX_CHROMOSOMES=sex,ANALYSIS_CHROMOSOMES=chromosomes)
    ref,sizes,views=reference(c)
    samples=read_table(samplesheet,('sample_id','matrix_path'))
    for s in samples:s['matrix_path']=resolve(s['matrix_path'],Path(samplesheet).parent)
    rows=[]
    for sample in samples:
        for res in resolutions:
            row={'sample_id':sample['sample_id'],'resolution':res,'status':'unavailable','uri':'','native_contacts':None,'analysis_bins':None,'finite_positive_weight_bins':None,'weight_converged':None,'reason':''}
            try:
                uri=matrix_uri(sample_matrix(c,sample,res),res);clr=cooler.Cooler(uri)
                stat=Path(uri.split('::')[0]).stat();row.update(file_size=stat.st_size,file_mtime_ns=stat.st_mtime_ns)
                bins=analysis_bins(sizes,views,res);ids=native_indices(clr,bins)
                require(clr.pixels()['count'][0:0].dtype.kind in 'iu','Raw counts are not integers')
                with clr.open('r') as group:
                    require('weight' in group['bins'],'Missing stored weights')
                    ds=group['bins/weight'];values=ds[:][ids];stats=dict(ds.attrs)
                finite=np.isfinite(values)
                require(not np.isinf(values).any() and (values[finite]>0).all() and finite.any(),'Invalid/absent usable weights')
                converged=bool(np.asarray(stats.get('converged',False)).all())
                row.update(uri=uri,native_contacts=int(clr.info['sum']),analysis_bins=len(ids),finite_positive_weight_bins=int(finite.sum()),weight_converged=converged,
                           status='metadata_usable' if converged else 'weight_convergence_unverified_or_failed',reason='Contact retention and biological reproducibility still require a pilot')
            except (ValueError,KeyError,OSError) as error:row['reason']=str(error)
            rows.append(row)
    destination.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(destination/'availability.tsv',sep='\t',index=False,na_rep='NA')
    atomic_json(destination/'audit.json',{'reference_id':c['REFERENCE_ID'],'reference_registry_sha256':file_hash(reference_registry),'scope':ref['effective_scope'],'samples':[s['sample_id'] for s in samples],
                                        'resolutions':resolutions,'rows':rows,'read_only':True,'eligibility_approved':False})
    print(destination/'availability.tsv')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--samples',required=True);p.add_argument('--reference-registry',required=True);p.add_argument('--out',required=True)
    p.add_argument('--resolutions',default='1000,2000,4000,10000,100000');p.add_argument('--matrix-registry',default='');p.add_argument('--chromosome-scope',default='autosomes');p.add_argument('--sex-chromosomes',default='');p.add_argument('--analysis-chromosomes',default='')
    a=p.parse_args();audit(a.samples,a.reference_registry,a.out,[int(x) for x in a.resolutions.split(',')],a.matrix_registry,a.chromosome_scope,a.sex_chromosomes,a.analysis_chromosomes)
