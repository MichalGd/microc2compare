#!/usr/bin/env python3
"""Verify current completion identities and required output hashes."""
import argparse
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from microc2compare.contracts import STAGES,enabled,read_json
from microc2compare.state import complete,scopes

def validate(c, until_stage='report'):
    failed=False
    for stage in STAGES[:STAGES.index(until_stage)+1]:
        if not enabled(c,stage):continue
        for scope in scopes(c,stage):
            valid=complete(c,stage,scope)
            print(f"{stage}/{scope}: {'valid' if valid else 'missing_or_stale'}")
            failed|=not valid
    if c['RUN_HICREP']:
        valid=complete(c,'hicrep','run');print(f"hicrep/run: {'valid' if valid else 'missing_or_stale'}");failed|=not valid
    return 1 if failed else 0


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config-json",required=True)
    p.add_argument("--until-stage",choices=STAGES,default="report")
    a=p.parse_args();sys.exit(validate(read_json(a.config_json),a.until_stage))
