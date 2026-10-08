#!/usr/bin/env python3
"""List shared Conda prefixes containing Rscript; never update a registry/env."""
import argparse
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--root", action="append", help="Shared prefix parent; repeat as needed")
a = p.parse_args()
found = set()
for root in a.root or ["/opt/conda_envs", "/opt/conda/envs", "/opt/miniconda/envs"]:
    parent = Path(root)
    if not parent.is_dir(): continue
    for prefix in [parent, *sorted(parent.iterdir())]:
        if (prefix/"conda-meta").is_dir() and (prefix/"bin/Rscript").is_file():
            found.add(str(prefix.resolve()))
for prefix in sorted(found): print(prefix)
if not found: print("No R Conda prefixes found in the supplied shared roots")
