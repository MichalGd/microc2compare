#!/usr/bin/env bash
declare -a M2C_CHILDREN=()

stop_children() {
  local pid
  for pid in "${M2C_CHILDREN[@]}"; do
    kill -TERM -- "-$pid" 2>/dev/null || kill -TERM "$pid" 2>/dev/null || true
  done
  for pid in "${M2C_CHILDREN[@]}"; do wait "$pid" 2>/dev/null || true; done
  M2C_CHILDREN=()
}

wait_batch() {
  local failed=0 pid
  for pid in "${M2C_CHILDREN[@]}"; do
    if wait "$pid"; then :; else failed=1; fi
  done
  M2C_CHILDREN=()
  [[ "$failed" == 0 ]]
}

batch_run() {
  local limit=$1; shift
  local callable=$1; shift
  local scope
  [[ "$limit" =~ ^[1-9][0-9]*$ ]] || die "Invalid concurrency: $limit"
  for scope in "$@"; do
    "$callable" "$scope" &
    M2C_CHILDREN+=("$!")
    if (( ${#M2C_CHILDREN[@]} >= limit )); then
      if wait_batch; then :; else stop_children; return 1; fi
    fi
  done
  wait_batch
}
