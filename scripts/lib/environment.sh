#!/usr/bin/env bash
env_run() {
  local kind=$1; shift
  local mode prefix key
  local -a measure=()
  if [[ -n "${M2C_TIME_OUTPUT:-}" && -x /usr/bin/time && "$(uname -s)" == Linux ]]; then
    measure=(/usr/bin/time -v -o "$M2C_TIME_OUTPUT")
  fi
  mode=$(config_get ENVIRONMENT_MODE)
  if [[ "$mode" == development ]]; then
    # Explicit local verification only; inventory still enters task signatures.
    if [[ "$1" == python ]]; then shift; "${measure[@]}" "$PYTHON_BOOTSTRAP" "$@";
    elif [[ "$1" == Rscript ]]; then shift; "${measure[@]}" "$R_SCRIPT" "$@";
    else "${measure[@]}" "$@"; fi
    return
  fi
  case "$kind" in core) key=CORE_ENV_PREFIX;; hicrep) key=HICREP_ENV_PREFIX;; diff_r) key=DIFF_R_ENV_PREFIX;; *) die "Unknown environment $kind";; esac
  prefix=$(config_get "$key")
  [[ -d "$prefix/conda-meta" ]] || die "Missing registered $kind prefix: $prefix; run explicit bootstrap setup"
  "${measure[@]}" "$CONDA_EXE" run --no-capture-output -p "$prefix" "$@"
}
