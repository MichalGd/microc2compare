#!/usr/bin/env bash
set -Eeuo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/lib/common.sh"
source "$M2C_ROOT/scripts/lib/environment.sh"
source "$M2C_ROOT/scripts/lib/parallel.sh"
source "$M2C_ROOT/scripts/lib/completion.sh"
set -m # Each asynchronous task has its own managed process group on Linux.
CONFIG_JSON= stage=
while (( $# )); do
  case "$1" in --config-json) CONFIG_JSON=$2; shift 2;; --stage) stage=$2; shift 2;; *) die "Unknown module option $1";; esac
done
[[ -n "$CONFIG_JSON" && -n "$stage" ]] || die "Module requires --config-json and --stage"
PYTHON_BOOTSTRAP=${PYTHON_BOOTSTRAP:-python3}
CONDA_EXE=$(config_get CONDA_EXE)
R_SCRIPT=$(config_get R_SCRIPT)
RUN_DIR=$(shell_path "$(config_get run_dir)")
trap stop_children EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
if ! state enabled --config "$CONFIG_JSON" --stage "$stage"; then
  state skipped --config "$CONFIG_JSON" --stage "$stage"
  exit 0
fi
if state restore-pools --config "$CONFIG_JSON" --stage "$stage"; then
  bash "$M2C_ROOT/scripts/run_stage.sh" --config-json "$CONFIG_JSON" --stage pool
fi
state prerequisites --config "$CONFIG_JSON" --stage "$stage"
limit=$(state concurrency --config "$CONFIG_JSON" --stage "$stage")
scope_text=$(state scopes --config "$CONFIG_JSON" --stage "$stage")
mapfile -t scopes <<< "$scope_text"
run_scope() {
  local scope=$1 kind=core
  if [[ "$stage" == chromosight || "$stage" == chromosight_discovery || "$stage" == loops && "$(config_get LOOP_CALLER)" == mustache ]]; then kind=loops; fi
  [[ "$stage" == stripes || "$stage" == stripe_consensus ]] && kind=stripes
  [[ "$stage" == browser_export && "$(config_get EXPORT_HIGLASS)" == true && -n "$(config_get BROWSER_ENV_PREFIX)" ]] && kind=browser
  [[ "$stage" == report ]] && kind=core
  local task_log="$RUN_DIR/logs/${stage}.${scope}.log"
  mkdir -p -- "$RUN_DIR/logs" "$RUN_DIR/benchmarks"
  export M2C_TIME_OUTPUT="$RUN_DIR/benchmarks/${stage}.${scope}.time.txt"
  log "$stage/$scope; concurrency=$limit; log=$task_log"
  env_run "$kind" python "$(python_path "$M2C_ROOT/scripts/scientific.py")" "$stage" --config "$CONFIG_JSON" --scope "$scope" > "$task_log" 2>&1
}
batch_run "$limit" run_scope "${scopes[@]}"
if [[ "$stage" == qc && "$(config_get RUN_HICREP)" == true ]]; then
  state concurrency --config "$CONFIG_JSON" --stage hicrep >/dev/null
  M2C_TIME_OUTPUT="$RUN_DIR/benchmarks/hicrep.run.time.txt" env_run hicrep python "$(python_path "$M2C_ROOT/scripts/scientific.py")" hicrep --config "$CONFIG_JSON" --scope run > "$RUN_DIR/logs/hicrep.run.log" 2>&1 &
  M2C_CHILDREN=("$!")
  wait_batch
fi
