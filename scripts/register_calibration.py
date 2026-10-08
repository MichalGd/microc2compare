#!/usr/bin/env python3
"""Register reviewed existing evidence after the actual preflight calibration gate.

No simulation rerun and no source/result edits; an explicit comparison/module
must match the current source, environment, design and recipe.
"""
import argparse
import sys
import tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from microc2compare.contracts import atomic_json, file_hash, read_json, require, validate_metadata
from microc2compare.science import check_calibration, environment_record


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("config-json", "results", "registry", "comparison-id", "module"): p.add_argument("--"+key,required=True)
    a=p.parse_args(); c=read_json(a.config_json); evidence=Path(a.results).resolve(); registry_path=Path(a.registry).resolve()
    require(evidence!=registry_path,"Registry and evidence must be different files")
    results=read_json(evidence)
    require(results.get("status")=="passed" and results.get("module")==a.module,"Evidence must pass and match the requested module")
    rows=[r for r in validate_metadata(c)["comparisons"] if r["comparison_id"]==a.comparison_id]
    require(len(rows)==1 and a.module in rows[0]["modules"].split(";"),"Comparison/module not present in configuration")
    comparison=dict(rows[0],modules=a.module)
    require(comparison["inference"]!="descriptive","Select an intended preliminary/formal design for registration")
    environments={"core":environment_record(c,"core"),"diff_r":environment_record(c,"diff_r")}
    if c["RUN_HICREP"]: environments["hicrep"]=environment_record(c,"hicrep")
    artifact={"calibration_key":results["calibration_key"],"status":"passed","results_path":str(evidence),"results_sha256":file_hash(evidence)}
    with tempfile.TemporaryDirectory(prefix="m2c-calibration-check-") as temporary:
        candidate=Path(temporary)/"registry.json"; atomic_json(candidate,{"schema_version":1,"artifacts":[artifact]})
        check_calibration(dict(c,CALIBRATION_REGISTRY=str(candidate),ALLOW_UNCALIBRATED_DEVELOPMENT=False),comparison,environments)
    registry=read_json(registry_path) if registry_path.exists() else {"schema_version":1,"artifacts":[]}
    if artifact not in registry["artifacts"]: registry["artifacts"].append(artifact)
    atomic_json(registry_path,registry)
    print(f"Registered matching reviewed evidence: {registry_path}")


if __name__=="__main__":
    try:main()
    except (ValueError,OSError) as error:print(f"ERROR: {error}",file=sys.stderr);sys.exit(2)
