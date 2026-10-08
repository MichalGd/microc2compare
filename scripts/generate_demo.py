#!/usr/bin/env python3
"""Generate original redistributable synthetic data; never modifies user maps."""
import argparse
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import cooler
import h5py
import numpy as np
import pandas as pd
from microc2compare.contracts import atomic_json,file_hash,write_table


def generate(output,development=False,raw_only=False,core_prefix=None):
    root=Path(output).resolve()
    if root.exists() and any(root.iterdir()):raise ValueError("Demo destination must be empty; choose a new output directory")
    root.mkdir(parents=True,exist_ok=True)
    sizes=root/"chrom.sizes";sizes.write_text("chrDemo\t6000000\n")
    views=root/"views.bed";views.write_text("chrDemo\t0\t6000000\tchrDemo\n")
    phase=root/"phasing.bedgraph"
    phase.write_text("".join(f"chrDemo\t{i*100000}\t{(i+1)*100000}\t{np.sin(i/7):.8f}\n" for i in range(60)))
    reference=root/"references.json"
    atomic_json(reference,{"reference_id":"demo","chromsizes":str(sizes),"chromsizes_sha256":file_hash(sizes),"views":str(views),"views_sha256":file_hash(views),"phasing_track":str(phase),"phasing_track_sha256":file_hash(phase)})
    bins=pd.DataFrame({"chrom":"chrDemo","start":np.arange(0,6000000,10000),"end":np.arange(10000,6000001,10000)})
    i,j=np.triu_indices(len(bins));distance=j-i
    compartment=np.sin(np.arange(len(bins))/70)
    means=60/(1+distance/10)*(1+0.4*np.sign(compartment[i]*compartment[j]))
    samples=[];membership=[]
    for index in range(6):
        sid=("C" if index<3 else "T")+str(index%3+1)
        condition="control" if index<3 else "treated"
        rng=np.random.default_rng(100+index)
        counts=rng.poisson(means*(1+index*.05)).astype(np.int64)
        focal=(i==100)&(j==180);counts[focal]+=5000 if index<3 else 10000
        pixels=pd.DataFrame({"bin1_id":i,"bin2_id":j,"count":counts})
        pixels=pixels.loc[pixels["count"]>0]
        matrix=root/(sid+".norm.mcool")
        cooler.create_cooler(str(matrix)+"::/resolutions/10000",bins,pixels,ordered=True,assembly="demo")
        cooler.coarsen_cooler(str(matrix)+"::/resolutions/10000",str(matrix)+"::/resolutions/100000",factor=10,chunksize=100000)
        balancing={}
        for resolution in (10000,100000):
            with h5py.File(matrix,"r+") as handle:
                group=handle[f"resolutions/{resolution}/bins"]
                group.create_dataset("weight",data=np.ones(len(group["start"])))
                group["weight"].attrs["converged"]=True
            balancing[str(resolution)]={"weight_column":"weight","convention":"multiplicative","converged":True,"recipe":{"kind":"known_synthetic_test_weights"}}
        provenance=root/(sid+".provenance.json")
        atomic_json(provenance,{"reference_id":"demo","count_units":"unique_contacts","filters":{"mapq":30,"min_distance_bp":0,"chromosome_policy":"chrDemo"},"balancing":balancing,"dedup_scope":"joint_library","library_ids":[sid+"_lib"],"source_ids":[sid+"_source"],"synthetic_fixture":True,"upstream_commit":"fixture_generator"})
        samples.append({"sample_id":sid,"project":"demo","condition":condition,"assay":"microc","biological_replicate":str(index%3+1),"subject_id":sid+"_culture","reference_id":"demo","matrix_path":str(matrix),"provenance_path":str(provenance),"batch":"batch"+str(index%3+1)})
        membership.append({"sample_id":sid,"library_id":sid+"_lib","source_id":sid+"_source","source_kind":"independent_library","dedup_group_id":sid+"_dedup","dedup_scope":"joint_library","source_checksum":file_hash(matrix)})
    write_table(root/"samples.csv",samples,list(samples[0]))
    write_table(root/"technical_membership.tsv",membership,list(membership[0]))
    comparison={"comparison_id":"treated_vs_control","case_condition":"treated","reference_condition":"control","design":"~condition","inference":"descriptive","modules":"loops" if raw_only else "loops;insulation"}
    write_table(root/"comparisons.tsv",[comparison],list(comparison))
    repo=Path(__file__).resolve().parents[1]
    config=root/"project.conf"
    lines=[f'source "{repo.as_posix()}/config/resources_demo.conf"',
           'PROJECT="demo"','ANALYSIS_SET_ID="synthetic_demo"','RUN_ID="run001"','REFERENCE_ID="demo"',
           'SAMPLESHEET="samples.csv"','TECHNICAL_MEMBERSHIP="technical_membership.tsv"','COMPARISONS="comparisons.tsv"',
           'REFERENCE_REGISTRY="references.json"','INDEPENDENT_FEATURES_REGISTRY=""','OUTDIR="results"','SCRATCH_ROOT="scratch"','RUN_POOL="false"','RUN_DEPTH_SENSITIVITY="false"']
    if core_prefix:lines.append('CORE_ENV_PREFIX="'+Path(core_prefix).as_posix()+'"')
    if development:lines += ['ENVIRONMENT_MODE="development"','PYTHON_BOOTSTRAP="'+Path(sys.executable).as_posix()+'"']
    if raw_only:lines += ['RUN_COMPARTMENTS="false"','RUN_INSULATION="false"','RUN_LOOPS="false"','RUN_AGGREGATES="false"','RUN_DIFFERENTIAL="false"']
    config.write_text("\n".join(lines)+"\n",encoding="utf-8")
    return config


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",required=True);p.add_argument("--development",action="store_true");p.add_argument("--raw-only",action="store_true");p.add_argument("--core-prefix")
    a=p.parse_args();print(generate(a.out,a.development,a.raw_only,a.core_prefix))
