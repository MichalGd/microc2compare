#!/usr/bin/env python3
"""Print the supported named analyses without requiring project data."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from microc2compare.analyses import CATALOGUE
for name, description in CATALOGUE.items():
    print(f"{name:24s} {description}")
