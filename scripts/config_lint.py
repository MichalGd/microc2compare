#!/usr/bin/env python3
"""Detect misspelled config assignments; trusted Bash is not a sandbox."""
import json
import re
import sys
from pathlib import Path

properties=json.loads((Path(__file__).resolve().parents[1]/"schemas/config.schema.json").read_text())["properties"]
text=Path(sys.argv[1]).read_text(encoding="utf-8-sig")
for line,content in enumerate(text.splitlines(),1):
    match=re.match(r"\s*(?:export\s+)?([A-Z][A-Z0-9_]*)=",content)
    if match and match[1] not in properties:
        raise SystemExit(f"Unknown configuration assignment {match[1]} at {sys.argv[1]}:{line}; check config_template.conf")
