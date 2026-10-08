#!/usr/bin/env python3
"""Check full/selected compartment agreement without requiring scientific packages."""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))
from microc2compare.acceptance import compare_compartments
from microc2compare.contracts import ContractError


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", required=True, help="Full-run sample .tsv; sibling eigenvalues/QC/bedGraph are checked")
    parser.add_argument("--selected", required=True, help="Selected-run sample .tsv")
    args = parser.parse_args()
    print(json.dumps(compare_compartments(args.full, args.selected), indent=2, allow_nan=False))


if __name__ == "__main__":
    try: main()
    except (ContractError, ValueError, KeyError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
