"""Optional fine-loop software is checked separately from the reused core."""
import sys
from pathlib import Path


def probe_command(c):
    script=Path(__file__).resolve().parents[2]/"scripts/caller_probe.py"
    options=[str(script),"--mustache" if c["RUN_LOOPS"] and c["LOOP_CALLER"]=="mustache" else "--no-mustache",
             "--chromosight" if c["RUN_CHROMOSIGHT"] or c.get("RUN_CHROMOSIGHT_DISCOVERY") else "--no-chromosight"]
    if c["ENVIRONMENT_MODE"]=="development":return [sys.executable,*options]
    return [c["CONDA_EXE"],"run","--no-capture-output","-p",c["LOOP_ENV_PREFIX"] or c["CORE_ENV_PREFIX"],"python",*options]
