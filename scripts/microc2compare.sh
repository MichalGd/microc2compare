#!/usr/bin/env bash
set -Eeuo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/lib/common.sh"
source "$M2C_ROOT/scripts/lib/environment.sh"
source "$M2C_ROOT/scripts/lib/parallel.sh"
source "$M2C_ROOT/scripts/lib/completion.sh"
set -m

usage() {
  cat <<'EOF'
Usage: bash scripts/microc2compare.sh -c PROJECT.conf [-s SAMPLES.csv]
       [--dry-run] [--resume] [--from-stage NAME] [--until-stage NAME]
       [--rerun-stage NAME]
       [--analysis NAME[,NAME...]] [--list-analyses]
       [--software-profile SOFTWARE.conf]
Reads existing biological matrices and stored weights; no sample rebalancing.
EOF
}
project_conf= sample_file= dry_run=false resume=false from_stage=preflight until_stage=report rerun_stage= analyses= list_analyses=false
while (( $# )); do
  case "$1" in
    -h|--help) usage; exit 0;;
    --list-analyses) list_analyses=true; shift;;
    -c|-s|--from-stage|--until-stage|--rerun-stage|--analysis|--software-profile)
      (( $# >= 2 )) || die "Missing value for $1"
      case "$1" in -c) project_conf=$2;; -s) sample_file=$2;; --from-stage) from_stage=$2;; --until-stage) until_stage=$2;; --rerun-stage) rerun_stage=$2;; --analysis) analyses=${analyses:+$analyses,}$2;; esac
      if [[ "$1" == --software-profile ]]; then
        [[ -n "$2" ]] || die "--software-profile requires a nonempty path"
        export M2C_SOFTWARE_PROFILE="$2"
      fi
      shift 2;;
    --dry-run) dry_run=true; shift;; --resume) resume=true; shift;;
    *) die "Unknown option: $1";;
  esac
done
if [[ "$list_analyses" == true ]]; then
  load_software_profile
  "${PYTHON_BOOTSTRAP:-python3}" "$(python_path "$M2C_ROOT/scripts/list_analyses.py")"
  exit 0
fi
[[ -n "$project_conf" ]] || { usage >&2; exit 2; }
[[ -z "$rerun_stage" || "$resume" == true ]] || die "--rerun-stage requires --resume"
export M2C_ANALYSES="$analyses"
load_config "$project_conf" "$sample_file"
M2C_LOCK_OWNED=false
finish() {
  local code=$?
  stop_children
  if [[ "$M2C_LOCK_OWNED" == true ]]; then
    if (( code != 0 )); then state failed --config "$CONFIG_JSON" --message "Runner exited with status $code" || true; fi
    rm -f -- "$RUN_DIR/.run-lock/owner.json"
    rmdir -- "$RUN_DIR/.run-lock" 2>/dev/null || true
  fi
  rm -f -- "$CONFIG_JSON"
}
trap finish EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
trap 'log "Command failed at line $LINENO"' ERR
selection=(--config "$CONFIG_JSON" --from-stage "$from_stage" --until-stage "$until_stage")
if [[ -n "$rerun_stage" ]]; then selection+=(--rerun-stage "$rerun_stage"); fi
state selection "${selection[@]}"
if [[ "$dry_run" == true ]]; then
  metadata plan --config "$CONFIG_JSON"
  state plan --config "$CONFIG_JSON" --from-stage "$from_stage" --until-stage "$until_stage"
  exit 0
fi
if [[ -d "$RUN_DIR" && "$resume" != true ]]; then die "Run directory exists; use --resume or a new RUN_ID"; fi
mkdir -p -- "$RUN_DIR"
if ! mkdir -- "$RUN_DIR/.run-lock" 2>/dev/null; then die "Run lock exists: $RUN_DIR/.run-lock; check owner liveness before explicit removal"; fi
M2C_LOCK_OWNED=true
state initialize --config "$CONFIG_JSON" --pid "$$" --hostname "$(hostname)"
if [[ -n "$rerun_stage" ]]; then state invalidate --config "$CONFIG_JSON" --stage "$rerun_stage"; fi
stage_text=$(state selection "${selection[@]}" --print)
mapfile -t stages <<< "$stage_text"
for stage in "${stages[@]}"; do
  log "Stage: $stage"
  bash "$M2C_ROOT/scripts/run_stage.sh" --config-json "$CONFIG_JSON" --stage "$stage" &
  M2C_CHILDREN=("$!")
  wait_batch
done
state summarize --config "$CONFIG_JSON"
log "Finished requested stages: $RUN_DIR"
