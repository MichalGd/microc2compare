#!/usr/bin/env bash
# Candidate-only acceptance: accepted pilot is read, never resumed or modified.
set -Eeuo pipefail
root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
site= baseline= destination=
while (( $# )); do
  case "$1" in
    --site) site=$2; shift 2;;
    --baseline-run) baseline=$2; shift 2;;
    --out) destination=$2; shift 2;;
    *) printf 'Usage: validate_performance_release.sh --site SOFTWARE.conf --baseline-run COMPLETED_PILOT --out NEW_DIRECTORY\n' >&2; exit 2;;
  esac
done
[[ -f "$site" && -f "$baseline/00_metadata/config.json" && -n "$destination" && ! -e "$destination" ]] || {
  printf 'Supply a software profile, completed pilot and new acceptance directory\n' >&2; exit 2;
}
umask 077
unset PYTHONHOME PYTHONPATH R_HOME R_LIBS R_LIBS_USER
source "$site"
mkdir -p "$destination"/{tmp,cache/matplotlib,cache/conda-pkgs,runtime-config,runtime-data}
export TMPDIR="$destination/tmp" XDG_CACHE_HOME="$destination/cache"
export XDG_CONFIG_HOME="$destination/runtime-config" XDG_DATA_HOME="$destination/runtime-data"
export MPLCONFIGDIR="$destination/cache/matplotlib" CONDA_PKGS_DIRS="$destination/cache/conda-pkgs"
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 CONDA_REGISTER_ENVS=false
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
"$CORE_ENV_PREFIX/bin/python" -I -B "$root/scripts/validate_tad_optimization.py" \
  --baseline-run "$baseline" --out "$destination/tad-replay" --workers 8
bash "$root/scripts/validate_fine_loop_demo.sh" --site "$site" --out "$destination/fine-demo"
printf 'Real-data TAD agreement and fine-loop demonstration passed\n' > "$destination/PASSED.txt"
printf 'Performance candidate acceptance passed. Audit: %s\n' "$destination"
