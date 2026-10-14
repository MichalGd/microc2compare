#!/usr/bin/env bash
declare -a M2C_CHILDREN=()
# Bash 5.1 added the PID result needed to remove exactly the job returned by
# wait -n. Older Bash versions use the same managed queue with short polling.
M2C_WAIT_N_P=0
if (( BASH_VERSINFO[0] > 5 || (BASH_VERSINFO[0] == 5 && BASH_VERSINFO[1] >= 1) )); then
  M2C_WAIT_N_P=1
fi

_m2c_remove_child() {
  local removed=$1 pid
  local -a remaining=()
  for pid in "${M2C_CHILDREN[@]}"; do
    [[ "$pid" == "$removed" ]] || remaining+=("$pid")
  done
  M2C_CHILDREN=("${remaining[@]}")
}

stop_children() {
  local pid
  for pid in "${M2C_CHILDREN[@]}"; do
    kill -TERM -- "-$pid" 2>/dev/null || kill -TERM "$pid" 2>/dev/null || true
  done
  for pid in "${M2C_CHILDREN[@]}"; do wait "$pid" 2>/dev/null || true; done
  M2C_CHILDREN=()
}

# Consume at least one completion, checking every already-finished child before
# returning a free slot. Bash can cache exit statuses for jobs which ended before
# wait -n: that builtin then returns 127 without setting its PID variable. Waiting
# for those managed PIDs explicitly recovers their real success/failure status.
_m2c_wait_one() {
  local pid completed_pid status completed=0
  while (( ${#M2C_CHILDREN[@]} )); do
    for pid in "${M2C_CHILDREN[@]}"; do
      if ! kill -0 "$pid" 2>/dev/null; then
        if wait "$pid"; then
          _m2c_remove_child "$pid"
          completed=1
        else
          _m2c_remove_child "$pid"
          return 1
        fi
      fi
    done
    (( completed == 0 )) || return 0
    if (( M2C_WAIT_N_P )); then
      completed_pid=
      if wait -n -p completed_pid "${M2C_CHILDREN[@]}" 2>/dev/null; then status=0; else status=$?; fi
      if [[ -n "${completed_pid:-}" ]]; then
        _m2c_remove_child "$completed_pid"
        (( status == 0 )) || return 1
        completed=1
      elif (( status != 127 )); then
        return 1
      else
        # All jobs may have ended between kill -0 and wait -n. The next queue
        # scan collects cached statuses; polling also covers older wait builtins.
        sleep 0.02
      fi
    else
      sleep 0.02
    fi
  done
  return 0
}

wait_batch() {
  while (( ${#M2C_CHILDREN[@]} )); do
    if _m2c_wait_one; then :; else stop_children; return 1; fi
  done
  return 0
}

batch_run() {
  local limit=$1; shift
  local callable=$1; shift
  local scope
  [[ "$limit" =~ ^[1-9][0-9]*$ ]] || die "Invalid concurrency: $limit"
  for scope in "$@"; do
    while (( ${#M2C_CHILDREN[@]} >= limit )); do
      if _m2c_wait_one; then :; else stop_children; return 1; fi
    done
    "$callable" "$scope" &
    M2C_CHILDREN+=("$!")
  done
  wait_batch
}
