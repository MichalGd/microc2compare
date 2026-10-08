#!/usr/bin/env bash
# Shared original helpers; sourced by the explicitly invoked entry/module scripts.
set -Eeuo pipefail
M2C_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)"
readonly M2C_ROOT
# Shared releases and scientific prefixes need no Python bytecode writes.
export PYTHONDONTWRITEBYTECODE=1

log() { printf '%s %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*" >&2; }
die() { log "ERROR: $*"; exit 2; }
python_path() { if command -v cygpath >/dev/null 2>&1; then cygpath -m "$1"; else printf '%s\n' "$1"; fi; }
shell_path() { if command -v cygpath >/dev/null 2>&1; then cygpath -u "$1"; else printf '%s\n' "$1"; fi; }
metadata() { "$PYTHON_BOOTSTRAP" "$(python_path "$M2C_ROOT/scripts/metadata.py")" "$@"; }

load_software_profile() {
  local profile=${M2C_SOFTWARE_PROFILE:-${M2C_DEFAULT_SOFTWARE_PROFILE:-}}
  local origin=launcher_default expected=${M2C_DEFAULT_SOFTWARE_PROFILE_SHA256:-} actual settings
  unset M2C_LOADED_SOFTWARE_PROFILE M2C_LOADED_SOFTWARE_PROFILE_SHA256 M2C_SOFTWARE_PROFILE_ORIGIN
  [[ -n "$profile" ]] || return 0
  if [[ -n "${M2C_SOFTWARE_PROFILE:-}" ]]; then origin=explicit_override; expected=; fi
  profile=$(shell_path "$profile")
  [[ -f "$profile" && -r "$profile" ]] || die "Software profile is not readable: $profile"
  profile="$(cd -- "$(dirname -- "$profile")" && pwd -P)/$(basename -- "$profile")"
  actual=$(sha256sum -- "$profile")
  actual=${actual%% *}
  [[ -z "$expected" || "$actual" == "$expected" ]] || die "Pinned software profile changed: $profile; publish a new versioned profile/launcher"
  # Profiles are trusted Bash, but only whitelisted software values cross this
  # child-shell boundary. Project paths, resources and analysis switches remain.
  # printf %q protects literal values; no profile stdout becomes executable text.
  settings=$(env -u BASH_ENV -u ENV bash --noprofile --norc -c '
    set -Eeuo pipefail
    source "$1" >&2
    for key in ENVIRONMENT_MODE CONDA_EXE PYTHON_BOOTSTRAP CORE_ENV_PREFIX HICREP_ENV_PREFIX DIFF_R_ENV_PREFIX CORE_LOCK HICREP_LOCK DIFF_R_LOCK R_SCRIPT; do
      if [[ -v "$key" ]]; then printf "export %s=%q\n" "$key" "${!key}"; fi
    done
  ' bash "$profile") || die "Software profile could not be loaded: $profile"
  eval "$settings"
  export M2C_LOADED_SOFTWARE_PROFILE="$(python_path "$profile")"
  export M2C_LOADED_SOFTWARE_PROFILE_SHA256="$actual" M2C_SOFTWARE_PROFILE_ORIGIN="$origin"
}

load_config() {
  local supplied=$1 sample_override=${2:-} invocation_dir=$PWD
  [[ -f "$supplied" ]] || die "Configuration not found: $supplied"
  CONFIG_PATH="$(cd -- "$(dirname -- "$supplied")" && pwd -P)/$(basename -- "$supplied")"
  CONFIG_DIR="$(dirname -- "$CONFIG_PATH")"
  set -a
  source "$M2C_ROOT/config/config_template.conf"
  pushd "$CONFIG_DIR" >/dev/null
  source "$CONFIG_PATH"
  popd >/dev/null
  load_software_profile
  if [[ -n "$sample_override" ]]; then
    if [[ "$sample_override" != /* && ! "$sample_override" =~ ^[A-Za-z]: ]]; then
      sample_override="$invocation_dir/$sample_override"
    fi
    SAMPLESHEET="$(python_path "$sample_override")"
  fi
  set +a
  "$PYTHON_BOOTSTRAP" "$(python_path "$M2C_ROOT/scripts/config_lint.py")" "$(python_path "$CONFIG_PATH")"
  export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  export NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
  CONFIG_JSON=$(mktemp)
  metadata normalize --config-dir "$(python_path "$CONFIG_DIR")" > "$CONFIG_JSON"
  CONFIG_JSON="$(python_path "$CONFIG_JSON")"
  RUN_DIR=$(shell_path "$(metadata get --config "$CONFIG_JSON" --key run_dir)")
}

config_get() { metadata get --config "$CONFIG_JSON" --key "$1"; }
