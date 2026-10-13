"""Exercise the real Bash scheduler using temporary files inside this workspace."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

import pytest


REPO = Path(__file__).resolve().parents[1]
BASH = (Path("C:/Program Files/Git/bin/bash.exe") if os.name == "nt"
        else Path(shutil.which("bash") or "/unavailable/bash"))
pytestmark = pytest.mark.skipif(not BASH.is_file(), reason="Bash unavailable")


@pytest.fixture
def workspace_tmp():
    parent = REPO / ".devtools" / "parallel-tests"
    parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=parent) as directory:
        yield Path(directory)


def run_shell(directory, text, *, poll=False):
    script = directory / "scheduler-test.sh"
    script.write_text("""#!/usr/bin/env bash
set -Eeuo pipefail
source "$1/scripts/lib/parallel.sh"
set -m
trap stop_children EXIT
die() { printf '%s\\n' "$*" >&2; exit 2; }
cd -- "$2"
""" + ("M2C_WAIT_N_P=0\n" if poll else "") + text, encoding="utf-8")
    environment = dict(os.environ, TMP=directory.as_posix(), TEMP=directory.as_posix(),
                       TMPDIR=directory.as_posix())
    started = time.monotonic()
    result = subprocess.run([str(BASH), "--noprofile", "--norc", script.as_posix(),
                             REPO.as_posix(), directory.as_posix()], cwd=REPO,
                            env=environment, capture_output=True, text=True, timeout=15)
    return result, time.monotonic() - started


@pytest.mark.parametrize("poll", [False, True], ids=["wait-n-p", "older-bash-poll"])
def test_free_slot_launches_third_while_first_still_running(workspace_tmp, poll):
    result, _ = run_shell(workspace_tmp, """
worker() {
  printf 'start %s\\n' "$1" >> events
  case "$1" in
    first)
      # A bounded barrier makes ordering independent of scheduler load.
      for ((i=0; i<200; i++)); do [[ -f third.started ]] && break; sleep 0.01; done;;
    second) sleep 0.1;;
    third) touch third.started; sleep 0.1;;
  esac
  printf 'end %s\\n' "$1" >> events
}
batch_run 2 worker first second third
[[ ${#M2C_CHILDREN[@]} == 0 ]]
""", poll=poll)
    assert result.returncode == 0, result.stderr
    events = (workspace_tmp / "events").read_text().splitlines()
    assert events.index("end second") < events.index("start third") < events.index("end first")


@pytest.mark.parametrize("poll", [False, True], ids=["wait-n-p", "older-bash-poll"])
def test_failure_stops_other_job_and_prevents_new_launch(workspace_tmp, poll):
    result, elapsed = run_shell(workspace_tmp, """
worker() {
  printf 'start %s\\n' "$1" >> events
  if [[ "$1" == slow ]]; then
    printf '%s\\n' "$BASHPID" > slow.pid
    exec sleep 10
  elif [[ "$1" == fail ]]; then
    sleep 0.1
    return 17
  fi
  printf 'end %s\\n' "$1" >> events
}
if batch_run 2 worker slow fail forbidden; then exit 8; fi
[[ ${#M2C_CHILDREN[@]} == 0 ]]
if kill -0 "$(cat slow.pid)" 2>/dev/null; then exit 9; fi
""", poll=poll)
    assert result.returncode == 0, result.stderr
    assert elapsed < 4
    assert "start forbidden" not in (workspace_tmp / "events").read_text()


@pytest.mark.parametrize("failure", [False, True], ids=["cached-success", "cached-failure"])
def test_wait_batch_recovers_status_when_every_child_already_ended(workspace_tmp, failure):
    result, _ = run_shell(workspace_tmp, """
worker() { return "$1"; }
worker 0 &
M2C_CHILDREN+=("$!")
worker """ + ("19" if failure else "0") + """ &
M2C_CHILDREN+=("$!")
sleep 0.2
if wait_batch; then status=0; else status=1; fi
[[ "$status" == """ + ("1" if failure else "0") + """ ]]
[[ ${#M2C_CHILDREN[@]} == 0 ]]
""")
    assert result.returncode == 0, result.stderr


def test_child_count_stays_bounded_across_many_fast_jobs(workspace_tmp):
    result, _ = run_shell(workspace_tmp, """
worker() {
  printf 'start %s\\n' "$1" >> events
  sleep 0.03
  printf 'end %s\\n' "$1" >> events
}
scopes=()
for ((i=0; i<30; i++)); do scopes+=("$i"); done
batch_run 3 worker "${scopes[@]}"
[[ ${#M2C_CHILDREN[@]} == 0 ]]
""")
    assert result.returncode == 0, result.stderr
    active, peak, finished = set(), 0, set()
    for line in (workspace_tmp / "events").read_text().splitlines():
        action, identity = line.split()
        if action == "start":
            assert identity not in active
            active.add(identity)
            peak = max(peak, len(active))
            assert len(active) <= 3
        else:
            active.remove(identity)
            finished.add(identity)
    assert peak == 3 and not active and len(finished) == 30


def test_wait_n_missing_pid_race_recovers_cached_failure(workspace_tmp):
    result, _ = run_shell(workspace_tmp, """
M2C_WAIT_N_P=1
simulate_race=1
wait() {
  if [[ "${1:-}" == -n && "$simulate_race" == 1 ]]; then
    simulate_race=0
    unset "$3"
    printf 'race\\n' > wait-race
    sleep 0.2
    return 127
  fi
  builtin wait "$@"
}
worker() { sleep 0.05; return "$1"; }
worker 0 &
M2C_CHILDREN+=("$!")
worker 19 &
M2C_CHILDREN+=("$!")
if wait_batch; then exit 8; fi
[[ ${#M2C_CHILDREN[@]} == 0 ]]
""")
    assert result.returncode == 0, result.stderr
    assert (workspace_tmp / "wait-race").read_text().strip() == "race"


@pytest.mark.skipif(os.name != "posix", reason="Linux process-group cleanup gate")
def test_failure_terminates_descendants_in_managed_group(workspace_tmp):
    result, elapsed = run_shell(workspace_tmp, """
worker() {
  if [[ "$1" == fail ]]; then sleep 0.1; return 23; fi
  trap 'kill -TERM "$descendant" 2>/dev/null || true; wait "$descendant" 2>/dev/null || true; exit 143' TERM
  sleep 10 &
  descendant=$!
  printf '%s\\n' "$descendant" > descendant.pid
  wait "$descendant"
}
if batch_run 2 worker slow fail; then exit 8; fi
if kill -0 "$(cat descendant.pid)" 2>/dev/null; then exit 9; fi
""")
    assert result.returncode == 0, result.stderr
    assert elapsed < 4
