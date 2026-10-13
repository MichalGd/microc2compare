#!/usr/bin/env bash
set -Eeuo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/lib/common.sh"
source "$M2C_ROOT/scripts/lib/environment.sh"
project_conf= action=
while (( $# )); do
  case "$1" in -c) project_conf=$2; shift 2;; --check-only|--create-missing) [[ -z "$action" ]] || die "Choose one bootstrap action"; action=$1; shift;; -h|--help) printf '%s\n' 'Usage: bootstrap_environments.sh -c PROJECT.conf (--check-only|--create-missing)'; exit 0;; *) die "Unknown option $1";; esac
done
[[ -n "$project_conf" && -n "$action" ]] || die "Supply configuration and exactly one action"
load_config "$project_conf"
trap 'rm -f -- "$CONFIG_JSON"' EXIT
required_text=$("$PYTHON_BOOTSTRAP" "$(python_path "$M2C_ROOT/scripts/environment_probe.py")" required --config "$CONFIG_JSON")
mapfile -t required <<< "$required_text"
for kind in "${required[@]}"; do
  case "$kind" in core) key=CORE_ENV_PREFIX; lockkey=CORE_LOCK; definition="$M2C_ROOT/environment.yml";; loops) key=LOOP_ENV_PREFIX; lockkey=LOOP_LOCK; definition="$M2C_ROOT/envs/fine-loops.yml";; browser) key=BROWSER_ENV_PREFIX; lockkey=BROWSER_LOCK; definition="$M2C_ROOT/envs/browser.yml";; stripes) key=STRIPE_ENV_PREFIX; lockkey=STRIPE_LOCK; definition="$M2C_ROOT/envs/stripes.yml";; hicrep) key=HICREP_ENV_PREFIX; lockkey=HICREP_LOCK; definition="$M2C_ROOT/envs/hicrep.yml";; diff_r) key=DIFF_R_ENV_PREFIX; lockkey=DIFF_R_LOCK; definition="$M2C_ROOT/envs/diff-r.yml";; esac
  prefix=$(config_get "$key")
  [[ "$kind" != stripes || -n "$prefix" ]] || die "Set a NEW STRIPE_ENV_PREFIX; the core environment will not be changed"
  if [[ "$kind" == loops && -z "$prefix" ]]; then
    [[ "$action" == --check-only ]] || die "Set a NEW LOOP_ENV_PREFIX before --create-missing; the core will not be changed"
    prefix=$(config_get CORE_ENV_PREFIX)
  fi
  if [[ "$(config_get ENVIRONMENT_MODE)" != development && ! -d "$prefix/conda-meta" ]]; then
    [[ "$action" == --create-missing ]] || die "Missing $kind prefix $prefix; set an owned path and run explicit --create-missing"
    [[ ! -e "$prefix" ]] || die "Refusing to install over an existing non-Conda directory: $prefix"
    lock=$(config_get "$lockkey")
    if [[ -n "$lock" ]]; then
      [[ -f "$lock" ]] || die "Configured lock missing: $lock"
      "$CONDA_EXE" create -y -p "$prefix" --file "$lock"
      if [[ "$kind" == hicrep || "$kind" == loops || "$kind" == stripes || "$kind" == browser ]]; then
        [[ -f "$lock.pip.txt" ]] || die "Exact $kind reconstruction needs a reviewed $lock.pip.txt with hashed pip requirements"
        env_run "$kind" python -m pip install --require-hashes --no-deps -r "$lock.pip.txt"
      fi
    else
      [[ "$(config_get ENVIRONMENT_MODE)" != locked ]] || die "Locked mode requires an exact tested lock; no unpinned fallback"
      "$CONDA_EXE" env create -y -p "$prefix" -f "$definition"
    fi
  fi
  "$PYTHON_BOOTSTRAP" "$(python_path "$M2C_ROOT/scripts/verify_environment_lock.py")" --config "$CONFIG_JSON" --kind "$kind"
  if [[ "$kind" == diff_r ]]; then
    env_run diff_r Rscript --vanilla -e 'for(p in c("edgeR","limma","statmod","jsonlite")) stopifnot(requireNamespace(p,quietly=TRUE)); sessionInfo()'
  else
    env_run "$kind" python "$(python_path "$M2C_ROOT/scripts/environment_probe.py")" "$kind" --config "$CONFIG_JSON"
  fi
done
