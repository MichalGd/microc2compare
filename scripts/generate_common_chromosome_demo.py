#!/usr/bin/env python3
"""Create NEW synthetic fixtures with X/Y and differing native chromosome scopes.

Only generated files below a previously empty destination are replaced. Real
matrices are never accepted as arguments. This is numerical/software validation,
not a biological simulation or statistical calibration.
"""
import argparse
import os
import sys
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parent))
from generate_demo import generate as original_demo
from microc2compare.contracts import atomic_json, file_hash, read_json, read_table, write_table


def generate(output, development=False, core_prefix=None, mixed=True):
    root=Path(output).resolve()
    configuration=original_demo(root,development=development,core_prefix=core_prefix)
    names=['chr1','chrX','chrY']
    sizes=root/'chrom.sizes'
    sizes.write_text(''.join(f'{chrom}\t2000000\n' for chrom in names)+'scaffold\t20000\n')
    views=root/'views.bed'
    views.write_text(''.join(f'{chrom}\t0\t2000000\t{chrom}\n' for chrom in names))
    phase=root/'phasing.bedgraph'
    phase.write_text(''.join(f'{chrom}\t{i*100000}\t{(i+1)*100000}\t{np.sin((i+k*20)/7):.8f}\n'
                             for k,chrom in enumerate(names) for i in range(20)))
    registry=read_json(root/'references.json')
    registry['chromosome_classes']={'autosomes':['chr1'],'sex':['chrX','chrY'],'other':['scaffold']}
    for key,path in [('chromsizes',sizes),('views',views),('phasing_track',phase)]:
        registry[key+'_sha256']=file_hash(path)
    atomic_json(root/'references.json',registry)
    selected=cooler.binnify(pd.Series({chrom:2000000 for chrom in names}),10000)
    for path in root.glob('*.norm.mcool'):
        sid=path.name.split('.')[0]
        raw=cooler.Cooler(str(path)+'::/resolutions/10000')
        pixels=raw.pixels()[:][['bin1_id','bin2_id','count']]
        heterogeneous=mixed and sid.startswith('C')
        native_names=['scaffold','chrY','chr1','chrX'] if heterogeneous else names
        native=cooler.binnify(pd.Series({chrom:20000 if chrom=='scaffold' else 2000000 for chrom in native_names}),10000)
        lookup=pd.MultiIndex.from_frame(native).get_indexer(pd.MultiIndex.from_frame(selected))
        a=lookup[pixels.bin1_id.to_numpy()];b=lookup[pixels.bin2_id.to_numpy()]
        pixels['bin1_id']=np.minimum(a,b);pixels['bin2_id']=np.maximum(a,b)
        if heterogeneous:
            pixels=pd.concat([pixels,pd.DataFrame({'bin1_id':[0,0,1],'bin2_id':[0,1,1],'count':[9000,5000,9000]})],ignore_index=True)
        pixels=pixels.sort_values(['bin1_id','bin2_id'])
        temporary=root/(sid+'.new.mcool')
        cooler.create_cooler(str(temporary)+'::/resolutions/10000',native,pixels,ordered=True,assembly='demo')
        cooler.coarsen_cooler(str(temporary)+'::/resolutions/10000',str(temporary)+'::/resolutions/100000',factor=10,chunksize=100000)
        with h5py.File(temporary,'r+') as handle:
            for resolution in (10000,100000):
                group=handle[f'resolutions/{resolution}/bins']
                ds=group.create_dataset('weight',data=np.ones(len(group['start'])))
                ds.attrs['converged']=True
        os.replace(temporary,path)
        prov=read_json(root/(sid+'.provenance.json'))
        prov['filters']['chromosome_policy']='all_with_scaffold' if heterogeneous else 'chr1_chrX_chrY'
        atomic_json(root/(sid+'.provenance.json'),prov)
    membership=read_table(root/'technical_membership.tsv')
    for row in membership:row['source_checksum']=file_hash(root/(row['sample_id']+'.norm.mcool'))
    write_table(root/'technical_membership.tsv',membership,list(membership[0]))
    configuration.write_text(configuration.read_text()+'\nANALYSIS_SET_ID="common_chromosomes"\nCHROMOSOME_SCOPE="autosomes_plus_sex"\nSEX_CHROMOSOMES="chrX,chrY"\n')
    return configuration


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',required=True);p.add_argument('--development',action='store_true')
    p.add_argument('--core-prefix');p.add_argument('--identical-scopes',action='store_true')
    args=p.parse_args()
    print(generate(args.out,args.development,args.core_prefix,not args.identical_scopes))
