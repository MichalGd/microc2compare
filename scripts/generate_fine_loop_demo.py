#!/usr/bin/env python3
"""Original six-replicate fine-grid fixture; not a biological calibration."""
import argparse
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import cooler
import h5py
import numpy as np
import pandas as pd
from generate_demo import generate as base_demo
from microc2compare.contracts import atomic_json,file_hash,read_json,read_table,write_table
from microc2compare.features import LOOP_COLUMNS


def generate(out,development=False):
    root=Path(out).resolve();config=base_demo(root,development=development)
    length=600000
    (root/'chrom.sizes').write_text(f'chrDemo\t{length}\n')
    (root/'views.bed').write_text(f'chrDemo\t0\t{length}\tchrDemo\n')
    registry=read_json(root/'references.json')
    for key,name in [('chromsizes','chrom.sizes'),('views','views.bed')]:registry[key+'_sha256']=file_hash(root/name)
    atomic_json(root/'references.json',registry)
    bins=cooler.binnify(pd.Series({'chrDemo':length}),2000)
    i,j=np.triu_indices(len(bins));distance=j-i
    compartment=np.sin(np.arange(len(bins))*2*np.pi/120)
    phase=root/'phasing.bedgraph'
    phase.write_text(''.join(f'chrDemo\t{start}\t{start+10000}\t{compartment[start//2000:(start+10000)//2000].mean():.8f}\n' for start in range(0,length,10000)))
    registry['phasing_track_sha256']=file_hash(phase)
    atomic_json(root/'references.json',registry)
    focal=np.exp(-((i-90)**2+(j-180)**2)/(2*1.2**2))
    for index,sid in enumerate(('C1','C2','C3','T1','T2','T3')):
        path=root/(sid+'.norm.mcool');path.unlink()
        means=80/(1+distance/8)*(1+0.6*compartment[i]*compartment[j])+focal*(3000 if index<3 else 6000)
        values=np.random.default_rng(1000+index).poisson(means).astype(np.int64)
        pixels=pd.DataFrame({'bin1_id':i,'bin2_id':j,'count':values});pixels=pixels.loc[pixels['count']>0]
        cooler.create_cooler(str(path)+'::/resolutions/2000',bins,pixels,ordered=True,assembly='demo')
        for res in (4000,10000,100000):cooler.coarsen_cooler(str(path)+'::/resolutions/2000',str(path)+f'::/resolutions/{res}',factor=res//2000,chunksize=100000)
        prov=read_json(root/(sid+'.provenance.json'))
        with h5py.File(path,'r+') as handle:
            for res in (2000,4000,10000,100000):
                group=handle[f'resolutions/{res}/bins'];ds=group.create_dataset('weight',data=np.ones(len(group['start'])));ds.attrs['converged']=True
                prov['balancing'][str(res)]={'weight_column':'weight','convention':'multiplicative','converged':True,'recipe':{'kind':'known_synthetic_test_weights'}}
        atomic_json(root/(sid+'.provenance.json'),prov)
    membership=read_table(root/'technical_membership.tsv')
    for row in membership:row['source_checksum']=file_hash(root/(row['sample_id']+'.norm.mcool'))
    write_table(root/'technical_membership.tsv',membership,list(membership[0]))
    # Analytical prespecified coordinates ensure scoring/quantification are exercised even if discovery is empty.
    candidates=root/'prespecified.bedpe.tsv'
    pd.DataFrame([['chrDemo',180000,182000,'chrDemo',360000,362000,'synthetic_focal']],columns=LOOP_COLUMNS).to_csv(candidates,sep='\t',index=False)
    atomic_json(root/'independent.json',{'feature_set_id':'analytic_fixture','reference_id':'demo','resolution':2000,'coordinate_convention':'0-based-half-open','path':str(candidates),'sha256':file_hash(candidates),'source':'original analytical synthetic generator; locations prespecified before drawing counts','discovery_samples':[],'independent_of_analysis_samples':True,'independence_evidence':'Coordinates are constants specified before synthetic count generation; no biological inference is claimed.'})
    gtf=root/'genes.gtf';gtf.write_text('chrDemo\tfixture\tgene\t180001\t210000\t.\t+\t.\tgene_id "geneA"; gene_name "DemoA";\nchrDemo\tfixture\texon\t180001\t182000\t.\t+\t.\tgene_id "geneA";\n')
    ccre=root/'ccre.bed';ccre.write_text('chrDemo\t360000\t362000\tDemoEnhancer\tdELS\n')
    atomic_json(root/'annotation.json',{'reference_id':'demo','gtf':str(gtf),'gtf_sha256':file_hash(gtf),'gtf_release':'original synthetic fixture','ccre':str(ccre),'ccre_sha256':file_hash(ccre),'ccre_source':'original synthetic intervals; not ENCODE data'})
    # The small fixture uses 10-kb compartments so its 600-kb view has >=50 bins.
    # Real-data compartments/saddles retain the 100-kb example setting.
    settings={'ANALYSIS_PROFILE':'practical_2_4kb','LOOP_CALLER':'mustache','LOOP_RESOLUTION':2000,'LOOP_DISCOVERY_RESOLUTIONS':'2000,4000','LOOP_QUANTIFICATION_RESOLUTIONS':'2000,4000','LOOP_MIN_DISTANCE_BP':5000,'LOOP_MAX_DISTANCE_BP':500000,'APA_FLANK_BP':100000,'RUN_CHROMOSIGHT':'true','RUN_ANNOTATION':'true','ANNOTATION_REGISTRY':'annotation.json','INDEPENDENT_FEATURES_REGISTRY':'independent.json','RUN_COMPARTMENTS':'true','COMPARTMENT_RESOLUTION':10000,'RUN_HICREP':'false','AGGREGATE_MODULES':'apa,ada,saddle','RUN_FINAL_REPORT':'true','ADA_RESOLUTION':10000,'TOTAL_CPUS':4,'MAX_MEMORY_MB':16000,'THREADS_LOOPS':1,'MEMORY_MB_LOOPS':8000}
    for kind in ('PREFLIGHT','QC','COMPARTMENTS','INSULATION','QUANTIFY','AGGREGATES','DIFFERENTIAL','REPORT'):settings.update({f'THREADS_{kind}':1,f'MEMORY_MB_{kind}':4000})
    with config.open('a') as handle:
        for key,value in settings.items():handle.write(f'\n{key}="{value}"')
        handle.write('\n')
    return config


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',required=True);p.add_argument('--development',action='store_true');a=p.parse_args();print(generate(a.out,a.development))
