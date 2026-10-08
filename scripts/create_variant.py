#!/usr/bin/env python3
"""Create a sibling project configuration with a new run ID and literal overrides."""
import argparse
import json
import shlex
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from microc2compare.contracts import IDENTIFIER, require


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config",required=True); p.add_argument("--out",required=True)
    p.add_argument("--run-id",required=True); p.add_argument("--set",action="append",default=[])
    a=p.parse_args(); source=Path(a.config).resolve(); target=Path(a.out).resolve()
    require(source.is_file(),"Source configuration missing")
    require(target.parent==source.parent,"Put the variant beside the original configuration to preserve relative paths")
    require(IDENTIFIER.fullmatch(a.run_id),"Supply a safe run ID")
    properties=json.loads((ROOT/"schemas/config.schema.json").read_text())["properties"]
    overrides={"RUN_ID":a.run_id}
    for assignment in a.set:
        key,separator,value=assignment.partition("=")
        require(separator and key in properties,"Use --set KNOWN_PARAMETER=literal_value")
        require(key!="RUN_ID","Use --run-id for the run identity")
        require("\n" not in value and "\r" not in value,"Parameter overrides must be single-line literals")
        overrides[key]=value
    text=source.read_text(encoding="utf-8-sig")+"\n# Variant: original settings retained; overrides below take precedence.\n"
    text+="\n".join(key+"="+shlex.quote(str(value)) for key,value in overrides.items())+"\n"
    with target.open("x",encoding="utf-8",newline="\n") as handle: handle.write(text)
    print(target)


if __name__=="__main__":
    try:main()
    except (ValueError,OSError) as error:print(f"ERROR: {error}",file=sys.stderr);sys.exit(2)
