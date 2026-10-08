import os
import signal
import subprocess
import time
from pathlib import Path

import pytest

REPO=Path(__file__).resolve().parents[1]


@pytest.mark.skipif(os.name!='posix',reason='Actual Linux process-group interruption gate')
def test_managed_interruption_reaps_child_group(tmp_path):
    script=tmp_path/'interrupt.sh';pidfile=tmp_path/'child.pid'
    script.write_text('''#!/usr/bin/env bash
set -Eeuo pipefail
source "$1/scripts/lib/common.sh"
source "$1/scripts/lib/parallel.sh"
set -m
trap stop_children EXIT
trap 'exit 143' TERM
worker(){ printf '%s\\n' "$BASHPID" > "$2"; exec sleep 30; }
worker ignored "$2" &
M2C_CHILDREN=("$!")
wait_batch
''')
    process=subprocess.Popen(['bash',str(script),str(REPO),str(pidfile)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,start_new_session=True)
    try:
        deadline=time.monotonic()+5
        while not pidfile.exists() and time.monotonic()<deadline:time.sleep(.02)
        assert pidfile.exists()
        child=int(pidfile.read_text())
        process.send_signal(signal.SIGTERM)
        process.communicate(timeout=5)
        assert process.returncode!=0
        with pytest.raises(ProcessLookupError):os.kill(child,0)
    finally:
        if process.poll() is None:
            os.killpg(process.pid,signal.SIGKILL);process.communicate(timeout=5)
