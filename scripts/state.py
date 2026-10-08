#!/usr/bin/env python3
import argparse
import json
import sys
from pathlib import Path
sys.stdout.reconfigure(newline="\n")

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from microc2compare.contracts import ContractError, STAGES, atomic_json, enabled, require
from microc2compare.state import (complete, initialize, invalidate, prerequisites, restore_required_pools,
                                  scopes, stage_concurrency)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("command", choices=("selection", "plan", "initialize", "invalidate", "enabled", "skipped", "prerequisites", "restore-pools", "concurrency", "scopes", "failed", "summarize"))
    p.add_argument("--config", required=True)
    p.add_argument("--stage")
    p.add_argument("--from-stage", default="preflight")
    p.add_argument("--until-stage", default="report")
    p.add_argument("--rerun-stage")
    p.add_argument("--print", action="store_true", dest="print_stages")
    p.add_argument("--pid", type=int)
    p.add_argument("--hostname")
    p.add_argument("--message")
    a = p.parse_args()
    c = json.loads(Path(a.config).read_text())
    if a.stage:
        require(a.stage in STAGES or a.stage == "hicrep", "Unknown stage")
    if a.command in ("selection", "plan"):
        require(a.from_stage in STAGES and a.until_stage in STAGES, "Unknown stage selection")
        require(STAGES.index(a.from_stage) <= STAGES.index(a.until_stage), "--from-stage follows --until-stage")
        require(not a.rerun_stage or a.rerun_stage in STAGES, "Unknown rerun stage")
        selected = STAGES[STAGES.index(a.from_stage):STAGES.index(a.until_stage)+1]
        if a.print_stages:
            print("\n".join(selected))
        if a.command == "plan":
            for s in selected:
                if enabled(c, s):
                    key = __import__('microc2compare.state', fromlist=['RESERVATIONS']).RESERVATIONS[s]
                    print(f"{s}: threads={c['THREADS_'+key]} memory_mb={c['MEMORY_MB_'+key]} concurrency={stage_concurrency(c,s)} prerequisites={','.join(prerequisites(c,s))}")
        return
    if a.command == "initialize": initialize(c, a.pid, a.hostname)
    elif a.command == "enabled": sys.exit(0 if enabled(c, a.stage) else 1)
    elif a.command == "concurrency": print(stage_concurrency(c, a.stage))
    elif a.command == "scopes": print("\n".join(scopes(c, a.stage)))
    elif a.command == "restore-pools": sys.exit(0 if restore_required_pools(c,a.stage) else 1)
    elif a.command == "invalidate": print("Invalidated: " + ",".join(invalidate(c, a.stage)))
    elif a.command == "prerequisites":
        for dep in prerequisites(c, a.stage):
            require(enabled(c, dep), f"{a.stage}: prerequisite {dep} is disabled")
            for ds in scopes(c, dep):
                require(complete(c, dep, ds), f"{a.stage}: {dep}/{ds} missing/stale; resume prerequisite")
    elif a.command == "skipped":
        atomic_json(Path(c["run_dir"]) / "00_metadata/status" / (a.stage+".json"), {"status": "skipped_config"})
    elif a.command in ("failed", "summarize"):
        records = [json.loads(x.read_text()) for x in (Path(c["run_dir"])/"completion").glob("*/*.done.json")]
        records=[r for r in records if enabled(c,r["stage"]) and (r["stage"]!="hicrep" or c["RUN_HICREP"])]
        for r in records:r["currently_valid"]=complete(c,r["stage"],r["scope"])
        atomic_json(Path(c["run_dir"])/"run_summary.json", {"status": "failed" if a.command == "failed" else "completed_requested_stages", "message": a.message, "tasks": records})


if __name__ == "__main__":
    try: main()
    except (ContractError, OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
