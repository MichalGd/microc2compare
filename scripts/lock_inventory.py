#!/usr/bin/env python3
"""Fingerprint a reviewed Conda package inventory; does not certify a lock."""
import argparse
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from microc2compare.contracts import atomic_json,digest,read_json
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--packages',required=True);p.add_argument('--out',required=True);a=p.parse_args()
atomic_json(a.out,{'package_digest':digest(read_json(a.packages)),'verification':'digest_only; attach actual Linux runtime verification evidence'})
