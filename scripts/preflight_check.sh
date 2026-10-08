#!/usr/bin/env bash
set -Eeuo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/lib/common.sh"
source "$M2C_ROOT/scripts/lib/environment.sh"
project_conf= sample_file=
while (( $# )); do
  case "$1" in -h|--help) printf '%s\n' 'Usage: preflight_check.sh -c PROJECT.conf [-s SAMPLES.csv]'; exit 0;; -c) project_conf=$2; shift 2;; -s) sample_file=$2; shift 2;; *) die "Unknown option $1";; esac
done
[[ -n "$project_conf" ]] || die "Supply -c PROJECT.conf"
load_config "$project_conf" "$sample_file"
trap 'rm -f -- "$CONFIG_JSON"' EXIT
env_run core python "$(python_path "$M2C_ROOT/scripts/scientific.py")" inspect --config "$CONFIG_JSON" --scope run
