#!/usr/bin/env python3
import argparse
import csv
import json
import os
import shutil
import tempfile
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from microc2compare.contracts import atomic_json, require
from microc2compare.state import complete, inventory, publish_completion, record_path


def validated_empty_browser_output(path):
    """Empty native tracks are valid only when the exporter inventory says so."""
    if path.suffix.lower() not in (".bed", ".bedpe", ".bedgraph"):
        return False
    for parent in path.parents:
        inventory_path = parent / "inventory.tsv"
        if inventory_path.is_file():
            with inventory_path.open(encoding="utf-8", newline="") as stream:
                for row in csv.DictReader(stream, delimiter="\t"):
                    if row.get("rows") == "0" and any(
                        row.get(key) and (parent / row[key]).resolve() == path.resolve()
                        for key in ("raw_path", "bedpe_path")
                    ):
                        return True
            return False
    return False


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("stage")
    parser.add_argument("--config",required=True)
    parser.add_argument("--scope",default="run")
    parser.add_argument("--derived-root")
    args=parser.parse_args()
    c=json.loads(Path(args.config).read_text())
    # Configure plotting before any scientific import. An explicit temporary
    # cache avoids writes to the user's home and Matplotlib's fallback atexit
    # cleanup racing with task scratch removal.
    runtime_cache=tempfile.TemporaryDirectory(prefix="microc2compare-runtime-")
    os.environ["MPLCONFIGDIR"]=runtime_cache.name
    from microc2compare.science import IMPLEMENTATIONS, cleanup_derived, inspect_inputs
    if args.stage=="inspect":
        print(json.dumps(inspect_inputs(c),indent=2));return
    if args.stage=="depth_hicrep":
        require(args.derived_root,"depth_hicrep requires a task-owned seed directory")
        target=Path(args.derived_root).resolve()
        require(target.is_relative_to((Path(c["run_dir"])/"01_qc/depth_sensitivity").resolve()),"Sensitivity path escapes task-owned root")
        from microc2compare.science import hicrep
        hicrep(c,args.scope,target);return
    if args.stage!="preflight" and complete(c,args.stage,args.scope):
        print(f"REUSED {args.stage}/{args.scope}");return
    # Remove the obsolete commit point before computing replacement products.
    done=record_path(c,args.stage,args.scope)
    if done.exists():done.unlink()
    started=time.monotonic()
    scratch=Path(c["SCRATCH_ROOT"])/c["PROJECT"]/c["ANALYSIS_SET_ID"]/c["RUN_ID"]/f"{args.stage}.{args.scope}"
    scratch.mkdir(parents=True,exist_ok=True)
    os.environ["TMPDIR"]=str(scratch)
    tempfile.tempdir=str(scratch)
    if args.stage=="report":
        from microc2compare.report import render
        outputs=render(c,args.scope)
    else:
        require(args.stage in IMPLEMENTATIONS,f"Unsupported stage {args.stage}")
        outputs=IMPLEMENTATIONS[args.stage](c,args.scope)
    require(outputs,f"{args.stage}: no validated product returned")
    for output in outputs:
        path=Path(output)
        require(path.is_file() and (path.stat().st_size > 0 or
                (args.stage == "browser_export" and validated_empty_browser_output(path))),
                f"Missing/empty required output {path}")
        if path.suffix==".json":json.loads(path.read_text())
    empty=False
    if args.stage in ("loops", "chromosight_discovery", "stripes"):
        import pandas as pd
        call_tables = [p for p in outputs if Path(p).suffix == ".tsv"]
        if args.stage == "stripes":
            call_tables = [p for p in call_tables if Path(p).name.startswith(args.scope + ".") and Path(p).name[len(args.scope)+1:-4].isdigit()]
        empty=bool(call_tables) and all(len(pd.read_csv(p,sep="\t"))==0 for p in call_tables)
    publish_completion(c,args.stage,args.scope,outputs,empty=empty,extra={"argv":sys.argv,"environment":inventory(),"elapsed_seconds":time.monotonic()-started,"exit_code":0,"normalization":"raw_counts" if args.stage in ("stripes", "stripe_consensus") else "export_of_saved_measurements" if args.stage == "browser_export" else "inherited_upstream_weights"})
    atomic_json(Path(c["run_dir"])/f"benchmarks/{args.stage}.{args.scope}.json",{"wall_seconds":time.monotonic()-started,"pid":os.getpid()})
    if args.stage=="report":cleanup_derived(c)
    if c["CLEAN_SCRATCH_ON_SUCCESS"]:
        require(scratch.resolve().is_relative_to(Path(c["SCRATCH_ROOT"]).resolve()),"Scratch cleanup escapes its registered root")
        shutil.rmtree(scratch)
    print(f"SUCCESS {args.stage}/{args.scope}: {len(outputs)} validated products")


if __name__=="__main__":
    try:main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
