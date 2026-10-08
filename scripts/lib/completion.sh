#!/usr/bin/env bash
# Python validates content hashes and dependency signatures; Bash only coordinates.
state() { "$PYTHON_BOOTSTRAP" "$(python_path "$M2C_ROOT/scripts/state.py")" "$@"; }
