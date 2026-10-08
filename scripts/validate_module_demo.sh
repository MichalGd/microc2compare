#!/usr/bin/env bash
# Linux acceptance in a new user-owned directory; imported fixture maps stay unchanged.
set -Eeuo pipefail
root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
site= destination= resource_config= resume=false
while (( $# )); do
  case "$1" in
    --site) site=$2; shift 2;;
    --out) destination=$2; shift 2;;
    --resource-config) resource_config=$2; shift 2;;
    --resume) resume=true; shift;;
    -h|--help) printf '%s\n' 'Usage: bash scripts/validate_module_demo.sh --site SOFTWARE.conf --out PRIVATE_DIRECTORY [--resource-config TRUSTED.conf] [--resume]'; exit 0;;
    *) printf 'Unknown option: %s\n' "$1" >&2; exit 2;;
  esac
done
[[ -f "$site" && -n "$destination" ]] || { printf '%s\n' 'Supply a readable software profile and private destination' >&2; exit 2; }
[[ "$resume" == true || ! -e "$destination" ]] || { printf '%s\n' 'Destination exists; use --resume for its unchanged acceptance fixture' >&2; exit 2; }
[[ -z "$resource_config" || -f "$resource_config" ]] || { printf '%s\n' 'Resource configuration is missing' >&2; exit 2; }
[[ "$resume" != true || -z "$resource_config" ]] || { printf '%s\n' 'Resume uses the existing resource configuration; do not supply a new resource preset' >&2; exit 2; }
umask 077
source "$site"
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
site="$(cd -- "$(dirname -- "$site")" && pwd -P)/$(basename -- "$site")"
resume_args=()
if [[ "$resume" == true ]]; then
  [[ -f "$destination/project.conf" && -s "$destination/input-sha256.txt" ]] || { printf '%s\n' 'Resume needs the original project.conf and input-sha256.txt' >&2; exit 2; }
  grep -Fxq -- "source \"$site\"" "$destination/project.conf" || { printf '%s\n' 'Existing fixture does not source this software profile; inspect configuration before resuming' >&2; exit 2; }
  sha256sum --check "$destination/input-sha256.txt"
  resume_args=(--resume)
else
  "$CONDA_EXE" run --no-capture-output -p "$CORE_ENV_PREFIX" python "$root/scripts/generate_demo.py" --out "$destination" --core-prefix "$CORE_ENV_PREFIX"
  if [[ -n "$resource_config" ]]; then
    resource_config="$(cd -- "$(dirname -- "$resource_config")" && pwd -P)/$(basename -- "$resource_config")"
    printf '\nsource "%s"\n' "$resource_config" >> "$destination/project.conf"
  fi
  printf '\nsource "%s"\n' "$site" >> "$destination/project.conf"
  sha256sum "$destination"/*.mcool > "$destination/input-sha256.txt"
fi
bash "$root/scripts/microc2compare.sh" -c "$destination/project.conf" --dry-run
bash "$root/scripts/microc2compare.sh" -c "$destination/project.conf" "${resume_args[@]}"
full="$destination/results/demo/synthetic_demo/run001"
"$PYTHON_BOOTSTRAP" "$root/scripts/validate_outputs.py" --config-json "$full/00_metadata/config.json"
bash "$root/scripts/microc2compare.sh" -c "$destination/project.conf" --resume
for analysis in qc hicrep compartments insulation tads loops quantify apa ada saddle insulation_comparison; do
  config="$destination/$analysis.conf"
  if [[ -e "$config" ]]; then
    [[ "$resume" == true && -f "$config" ]] || { printf 'Refusing to overwrite %s\n' "$config" >&2; exit 2; }
  else
    cp -- "$destination/project.conf" "$config"
    printf '\nRUN_ID="%s-only"\n' "$analysis" >> "$config"
  fi
  bash "$root/scripts/microc2compare.sh" -c "$config" --analysis "$analysis" "${resume_args[@]}"
  run="$destination/results/demo/synthetic_demo/$analysis-only"
  "$PYTHON_BOOTSTRAP" "$root/scripts/validate_outputs.py" --config-json "$run/00_metadata/config.json"
  if [[ "$analysis" == insulation || "$analysis" == loops || "$analysis" == compartments ]]; then
    for sample in C1 C2 C3 T1 T2 T3; do
      case "$analysis" in
        insulation) product="03_features/insulation/$sample.tsv";;
        loops) product="03_features/loops/$sample.10000.bedpe.tsv";;
        compartments) product="03_features/compartments/$sample.tsv";;
      esac
      if [[ "$analysis" == compartments ]]; then
        "$PYTHON_BOOTSTRAP" "$root/scripts/compare_compartments.py" --full "$full/$product" --selected "$run/$product"
      else
        cmp -- "$full/$product" "$run/$product"
      fi
    done
  fi
done
if [[ -e "$destination/tad-variant.conf" ]]; then
  [[ "$resume" == true && -f "$destination/tad-variant.conf" ]] || { printf '%s\n' 'Refusing to overwrite TAD variant configuration' >&2; exit 2; }
else
  cp -- "$destination/project.conf" "$destination/tad-variant.conf"
  printf '\nRUN_ID="tad-threshold-variant"\nINSULATION_BOUNDARY_THRESHOLD="0.5"\n' >> "$destination/tad-variant.conf"
fi
bash "$root/scripts/microc2compare.sh" -c "$destination/tad-variant.conf" --analysis tads "${resume_args[@]}"
"$PYTHON_BOOTSTRAP" "$root/scripts/validate_outputs.py" --config-json "$destination/results/demo/synthetic_demo/tad-threshold-variant/00_metadata/config.json"
sha256sum --check "$destination/input-sha256.txt"
[[ -z "$(find "$destination/results" -type f \( -name '*.cool' -o -name '*.mcool' \) -print -quit)" ]]
printf 'Full/selected/variant acceptance passed. Report: %s/08_reports/index.html\n' "$full"
