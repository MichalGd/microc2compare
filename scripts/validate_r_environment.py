#!/usr/bin/env python3
"""Read-only R package probe and numerical adapter acceptance in a new audit dir.

Run using the checked core Python. This never installs packages or enables
formal inference; passing numerical checks do not replace calibration.
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from microc2compare.contracts import atomic_json, file_hash, require


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prefix", required=True)
    p.add_argument("--conda", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    prefix = Path(a.prefix).resolve(); out = Path(a.out).resolve()
    require((prefix/"conda-meta").is_dir(), "Supply an existing Conda prefix; no creation occurs here")
    require(not out.exists(), "Choose a new private audit output directory")
    require(not out.is_relative_to(prefix), "Audit outputs must be outside the tested environment")
    out.mkdir(parents=True, mode=0o700)
    env = dict(os.environ, R_ENVIRON_USER="/dev/null", R_PROFILE_USER="/dev/null")
    for key in ("R_HOME", "R_LIBS", "R_LIBS_USER"): env.pop(key, None)
    command = [a.conda, "run", "--no-capture-output", "-p", str(prefix), "Rscript", "--vanilla"]
    def execute(name, args):
        result = subprocess.run(command + list(map(str, args)), capture_output=True, text=True, env=env)
        (out/(name+".log")).write_text(result.stdout+"\n"+result.stderr, encoding="utf-8")
        require(result.returncode == 0, f"{name} failed; inspect {out/(name+'.log')}")
        return result.stdout
    evidence = {"prefix": str(prefix), "status": "failed", "calibration_status": "not_run",
                "adapters": {n: file_hash(ROOT/"scripts"/n) for n in ("differential_loops.R", "differential_loop_uncertainty.R", "differential_insulation.R")}}
    try:
        runtime = execute("packages", ["-e", 'for(p in c("edgeR","limma","statmod","jsonlite")) stopifnot(requireNamespace(p,quietly=TRUE)); sessionInfo()'])
        evidence["session_info"] = runtime
        rng = np.random.default_rng(97); ids = ["C1", "C2", "C3", "T1", "T2", "T3"]; n = 1200
        params = {"design": {"sample_ids": ids, "matrix": [[1,0]]*3+[[1,1]]*3, "contrast": [0,1]},
                  "validation_status": "uncalibrated", "inference": "formal", "fdr": .05,
                  "min_total": 20, "min_mean_at_min_exposure": 5, "min_features": 1000}
        atomic_json(out/"parameters.json", params)
        means = np.full((n,6), 80.); means[:60,3:] *= 3
        counts = pd.DataFrame(rng.negative_binomial(20,20/(20+means)), columns=ids)
        counts.insert(0, "feature_id", [f"f{i}" for i in range(n)])
        counts.to_csv(out/"counts.tsv", sep="\t", index=False)
        ex = pd.DataFrame({"sample_id": ids, "eligible_cis_exposure": [1000000]*6})
        ex.to_csv(out/"exposures.tsv", sep="\t", index=False)
        execute("loops", [ROOT/"scripts/differential_loops.R", out/"counts.tsv", out/"exposures.tsv", out/"parameters.json", out/"loops.tsv"])
        loops = pd.read_csv(out/"loops.tsv", sep="\t").set_index("feature_id")
        affected = [f"f{i}" for i in range(60)]
        require(loops.loc[affected, "logFC"].median() > 1, "Loop effect direction/magnitude failed")
        require(loops.loc[[f"f{i}" for i in range(60,n)], "logFC"].abs().median() < .2, "Loop null effects failed")
        require(loops.FDR.between(0,1).all() and loops.SE_log2_approx_QL.gt(0).all(), "Invalid loop uncertainty/FDR")
        ex.eligible_cis_exposure *= 7; ex.to_csv(out/"exposures-scaled.tsv", sep="\t", index=False)
        execute("loops-scaled", [ROOT/"scripts/differential_loops.R", out/"counts.tsv", out/"exposures-scaled.tsv", out/"parameters.json", out/"loops-scaled.tsv"])
        scaled = pd.read_csv(out/"loops-scaled.tsv", sep="\t").set_index("feature_id")
        np.testing.assert_allclose(loops.logFC, scaled.loc[loops.index,"logFC"], atol=.005)
        scores = pd.DataFrame(rng.normal(0,.1,(n,6)), columns=ids); scores.loc[:59,ids[3:]] += 1
        scores.loc[60,"C1"] = np.nan; scores.insert(0,"feature_id",[f"f{i}" for i in range(n)])
        scores.to_csv(out/"scores.tsv", sep="\t", index=False)
        execute("insulation", [ROOT/"scripts/differential_insulation.R", out/"scores.tsv", out/"parameters.json", out/"insulation.tsv"])
        ins = pd.read_csv(out/"insulation.tsv", sep="\t").set_index("feature_id")
        require(ins.loc[affected,"case_minus_reference_log2_insulation"].median() > .8, "Insulation effect direction failed")
        require(pd.isna(ins.loc["f60","PValue"]) and pd.isna(ins.loc["f60","FDR"]), "Insulation missing value was incorrectly tested")
        def bh(values):
            values = np.asarray(values); order = np.argsort(values); ranked = values[order]*len(values)/np.arange(1,len(values)+1)
            adjusted = np.minimum(1,np.minimum.accumulate(ranked[::-1])[::-1]); result = np.empty(len(values)); result[order] = adjusted
            return result
        for frame in (loops, ins):
            good = frame.PValue.notna(); np.testing.assert_allclose(frame.loc[good,"FDR"], bh(frame.loc[good,"PValue"]), atol=1e-10)
        evidence["status"] = "passed_numerical_checks_only"
        print("R package and numerical checks passed; statistical calibration remains required")
    finally:
        atomic_json(out/"acceptance.json", evidence)


if __name__ == "__main__":
    try: main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr); sys.exit(1)
