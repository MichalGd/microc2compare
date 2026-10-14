#!/usr/bin/env bash
set -Eeuo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/lib/common.sh"
source "$M2C_ROOT/scripts/lib/environment.sh"
source "$M2C_ROOT/scripts/lib/parallel.sh"
source "$M2C_ROOT/scripts/lib/completion.sh"
set -m

usage() {
  cat <<'EOF'
Usage: microc2compare -c PROJECT.conf [-s SAMPLES.csv]
       [--dry-run] [--resume] [--from-stage NAME] [--until-stage NAME]
       [--rerun-stage NAME]
       [--analysis NAME[,NAME...]] [--list-analyses]
       [--software-profile SOFTWARE.conf]
       --init PROJECT.conf [--preset existing|practical_2_4kb|hsieh_8grid|short_range_4kb] [--reference ASSEMBLY]
       -c PROJECT.conf --prepare-inputs|--register-resources|--audit-resolutions
       -c PROJECT.conf --prepare-resolutions|--validate|--demo [--out NEW_DIRECTORY]
       -c PROJECT.conf --prepare-from-pairs|--input-capabilities|--export-browser [--out NEW_DIRECTORY]
       -c PROJECT.conf --variant NEW_RUN_ID [--set PARAMETER=value ...]
Reads existing biological matrices and stored weights; no sample rebalancing.
Preparation commands are explicit, separate actions; --prepare-resolutions
balances only newly created exact coarsenings, never imported matrices.
EOF
}
project_conf= sample_file= dry_run=false resume=false from_stage=preflight until_stage=report rerun_stage= analyses= list_analyses=false
action= action_out= preset=practical_2_4kb init_conf= preset_supplied=false
reference_id=mm39 reference_supplied=false
variant_run=; variant_sets=()
set_action() { [[ -z "$action" ]] || die "Choose one preparation/checking action"; action=$1; }
while (( $# )); do
  case "$1" in
    -h|--help) usage; exit 0;;
    --list-analyses) list_analyses=true; shift;;
    --prepare-inputs|--register-resources|--audit-resolutions|--prepare-resolutions|--prepare-from-pairs|--input-capabilities|--export-browser|--validate|--demo)
      set_action "${1#--}"; shift;;
    --variant|--set)
      (( $# >= 2 )) && [[ -n "$2" ]] || die "Missing value for $1"
      case "$1" in --variant) set_action variant; variant_run=$2;; --set) variant_sets+=("$2");; esac
      shift 2;;
    --init|--preset|--out|--reference)
      (( $# >= 2 )) && [[ -n "$2" ]] || die "Missing value for $1"
      case "$1" in --init) set_action init; init_conf=$2;; --preset) preset=$2; preset_supplied=true;; --out) action_out=$2;; --reference) reference_id=$2; reference_supplied=true;; esac
      shift 2;;
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
if [[ -n "$action" ]]; then
  [[ "$dry_run" == false && "$resume" == false && -z "$analyses" && "$list_analyses" == false && -z "$rerun_stage" && "$from_stage" == preflight && "$until_stage" == report && -z "$sample_file" ]] || die "Preparation/checking actions cannot be combined with analysis execution options"
fi
[[ -z "$action_out" || ( -n "$action" && "$action" != init && "$action" != validate && "$action" != variant ) ]] || die "--out applies to preparation/audit actions"
[[ "$preset_supplied" == false || "$action" == init ]] || die "--preset applies only to --init"
[[ "$reference_supplied" == false || "$action" == init ]] || die "--reference applies only to --init"
(( ${#variant_sets[@]} == 0 )) || [[ "$action" == variant ]] || die "--set applies only to --variant"
if [[ "$action" == init ]]; then
  [[ -z "$project_conf" ]] || die "Use --init PROJECT.conf without -c"
  load_software_profile
  "${PYTHON_BOOTSTRAP:-python3}" "$(python_path "$M2C_ROOT/scripts/project_tools.py")" --action init --project-config "$(python_path "$init_conf")" --preset "$preset" --reference "$reference_id"
  exit 0
fi
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
if [[ -n "$action" ]]; then
  action_args=(--action "$action" --config-json "$CONFIG_JSON" --project-config "$(python_path "$CONFIG_PATH")")
  if [[ -n "$action_out" ]]; then action_args+=(--out "$(python_path "$action_out")"); fi
  if [[ "$action" == variant ]]; then
    action_args+=(--run-id "$variant_run")
    for setting in "${variant_sets[@]}"; do action_args+=(--set "$setting"); done
  fi
  action_kind=core
  [[ "$action" == export-browser && "$(config_get EXPORT_HIGLASS)" == true && -n "$(config_get BROWSER_ENV_PREFIX)" ]] && action_kind=browser
  env_run "$action_kind" python "$(python_path "$M2C_ROOT/scripts/project_tools.py")" "${action_args[@]}"
  exit 0
fi
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
