#!/usr/bin/env python3
"""Run the shipped R models on frozen independent simulation scenarios."""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np
import pandas as pd
from scipy.signal import lfilter
from scipy.stats import beta
from microc2compare.contracts import atomic_json,file_hash,read_json,require
from microc2compare.science import calibration_design_domain,calibration_key,design_matrix,environment_record,run_r,save
from microc2compare.state import source_digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config-json",required=True,help="Private normalized config produced by an installed workflow")
    p.add_argument("--module",required=True,choices=("loops","insulation"))
    p.add_argument("--out",required=True)
    p.add_argument("--datasets",type=int,default=200)
    p.add_argument("--seed",type=int,default=20261007)
    p.add_argument("--design",choices=("unpaired","paired"),default="unpaired")
    p.add_argument("--units-per-condition",type=int,default=3)
    p.add_argument("--register",help="Explicitly write a calibration registry entry after all criteria pass")
    a=p.parse_args();require(a.datasets>=200,"Primary calibration requires at least 200 independent datasets")
    require(a.units_per_condition>=2,"Calibration requires at least two independent units per condition")
    c=read_json(a.config_json);out=Path(a.out)
    require(not out.exists(),"Calibration output exists; choose a new directory to preserve evidence")
    out.mkdir(parents=True)
    require(c["FDR_THRESHOLD"]==.05,"The prespecified v0.1 primary protocol fixes nominal FDR at 0.05; revise/review the protocol before testing other thresholds")
    environments={"core":environment_record(c,"core"),"diff_r":environment_record(c,"diff_r")}
    if c["RUN_HICREP"]:environments["hicrep"]=environment_record(c,"hicrep")
    scenarios=["equal_depth_null","unequal_depth_null","mixture_10pct","broad_30pct","sparse_null"]
    summaries=[]
    units=a.units_per_condition;n_samples=2*units
    comparison={"comparison_id":"simulation","reference_condition":"control","case_condition":"treated","design":"~subject_id+condition" if a.design=="paired" else "~condition","inference":"formal" if units>=3 else "preliminary","modules":a.module}
    samples=[{"sample_id":"S"+str(i),"condition":"control" if i<units else "treated","subject_id":"donor"+str(i%units if a.design=="paired" else i)} for i in range(n_samples)]
    if a.design=="paired":scenarios.append("paired_null")
    design=design_matrix(samples,comparison,c);ids=design["sample_ids"]
    domain=calibration_design_domain(c,comparison,samples)
    key=calibration_key(c,a.module,comparison,environments,design_domain=domain)
    params={"design":design,"validation_status":"simulation_calibration","inference":"formal","fdr":c["FDR_THRESHOLD"],"min_total":c["LOOP_MIN_TOTAL_COUNT"],"min_mean_at_min_exposure":c["LOOP_MIN_MEAN_COUNT_AT_MIN_EXPOSURE"],"min_features":c["LOOP_MIN_FEATURES_FORMAL"]}
    atomic_json(out/"parameters.json",params)
    for scenario_index,scenario in enumerate(scenarios):
        fractions=[];powers=[]
        for iteration in range(a.datasets):
            rng=np.random.default_rng(np.random.SeedSequence([a.seed,scenario_index,iteration]))
            features=max(2500,2*c["LOOP_MIN_FEATURES_FORMAL"])
            truth=np.zeros(features,dtype=bool)
            if scenario=="mixture_10pct":truth[:features//10]=True
            if scenario=="broad_30pct":truth[:3*features//10]=True
            if a.module=="loops":
                depth=np.ones(n_samples)
                if scenario=="unequal_depth_null":depth=np.linspace(.6,2,n_samples)
                base=rng.lognormal(np.log(40),.5,features)
                if scenario=="sparse_null":base[:features//5]=.1
                means=base[:,None]*depth[None,:]
                means[truth,units:]*=2
                if a.design=="paired":means*=np.tile(rng.lognormal(0,.2,(features,units)),(1,2))
                dispersion=.1
                y=rng.negative_binomial(1/dispersion,1/(1+means*dispersion))
                # Prespecified global eligible totals define relative-allocation truth.
                # Background contacts fill the remaining eligible domain, so a broad
                # effect does not silently change every null feature's denominator.
                exposures=np.rint(10000000*depth).astype(np.int64)
                require((exposures-y.sum(axis=0)>0).all(),"Candidate contacts exhaust the eligible background domain")
                frame=pd.DataFrame(y,columns=ids);frame.insert(0,"feature_id",["loop_"+str(i) for i in range(features)])
                save(frame,out/"input.tsv");save(pd.DataFrame({"sample_id":ids,"eligible_cis_exposure":exposures}),out/"exposures.tsv")
                run_r(c,"differential_loops.R",[out/"input.tsv",out/"exposures.tsv",out/"parameters.json",out/"result.tsv"])
            else:
                mean=rng.normal(0,.5,features);noise=lfilter([1],[1,-.6],rng.normal(0,.25,(features,n_samples)),axis=0)
                scores=mean[:,None]+noise;scores[truth,units:]+=.5
                if a.design=="paired":
                    subject_effect=rng.normal(0,.4,(features,units))
                    scores+=np.tile(subject_effect,(1,2))
                frame=pd.DataFrame(scores,columns=ids);frame.insert(0,"feature_id",["bin_"+str(i) for i in range(features)])
                save(frame,out/"input.tsv");run_r(c,"differential_insulation.R",[out/"input.tsv",out/"parameters.json",out/"result.tsv"])
            result=pd.read_csv(out/"result.tsv",sep="\t")
            selected=set(result.loc[result.FDR<=c["FDR_THRESHOLD"],"feature_id"])
            true_ids=set(frame.loc[truth,"feature_id"])
            fractions.append(len(selected-true_ids)/max(1,len(selected)))
            powers.append(len(selected&true_ids)/max(1,len(true_ids)))
        bootstrap=np.random.default_rng(a.seed+scenario_index).choice(fractions,size=(5000,len(fractions)),replace=True).mean(axis=1)
        upper=max(float(np.quantile(bootstrap,.975)),float(np.mean(fractions))+1.96*np.std(fractions,ddof=1)/np.sqrt(len(fractions)))
        if not any(fractions):upper=float(beta.ppf(.975,1,len(fractions)))
        summaries.append({"scenario":scenario,"n_datasets":a.datasets,"mean_fdp":float(np.mean(fractions)),"fdp_upper_95":upper,"mean_power":float(np.mean(powers))})
        print(json.dumps(summaries[-1]),flush=True)
    results={"calibration_key":key,"source_digest":source_digest("differential"),"environments":environments,"module":a.module,"seed":a.seed,"design_domain":domain,"primary_scenarios":summaries,"status":"passed" if all(x["fdp_upper_95"]<=.075 for x in summaries) else "failed","scope":"simulated_"+a.design+"_design_only; biological negative controls pending"}
    path=out/"calibration_results.json";atomic_json(path,results)
    if a.register:
        require(results["status"]=="passed","Failed suite cannot be registered as passing")
        registry_path=Path(a.register);registry=read_json(registry_path) if registry_path.exists() else {"schema_version":1,"artifacts":[]}
        registry["artifacts"].append({"calibration_key":key,"status":"passed","results_path":str(path.resolve()),"results_sha256":file_hash(path)})
        atomic_json(registry_path,registry)
    require(results["status"]=="passed","Statistical calibration failed; do not enable production inference")


if __name__=="__main__":main()
