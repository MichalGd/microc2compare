#!/usr/bin/env bash
# Linux acceptance using actual optional caller software; no biological calibration claim.
set -Eeuo pipefail
root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
site= destination= resource= resume=false
while (( $# )); do
  case "$1" in
    --site) site=$2; shift 2;; --out) destination=$2; shift 2;; --resource-config) resource=$2; shift 2;; --resume) resume=true; shift;;
    *) printf 'Usage: validate_fine_loop_demo.sh --site SOFTWARE.conf --out NEW_PRIVATE_DIRECTORY [--resume]\n' >&2; exit 2;;
  esac
done
[[ -f "$site" && -n "$destination" ]] || { printf 'Missing software profile/destination\n' >&2; exit 2; }
[[ -z "$resource" || -f "$resource" && "$resume" == false ]] || { printf 'Resource config is only accepted for a new fixture\n' >&2; exit 2; }
umask 077
unset PYTHONHOME PYTHONPATH R_HOME R_LIBS R_LIBS_USER
source "$site"
[[ -n "${LOOP_ENV_PREFIX:-}" ]] || { printf 'Register LOOP_ENV_PREFIX first\n' >&2; exit 2; }
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
if [[ "$resume" == false ]]; then
  [[ ! -e "$destination" ]] || { printf 'Destination exists; choose a new acceptance directory\n' >&2; exit 2; }
  "$CORE_ENV_PREFIX/bin/python" -I -B "$root/scripts/generate_fine_loop_demo.py" --out "$destination"
  printf '\nsource "%s"\n' "$site" >> "$destination/project.conf"
  if [[ -n "$resource" ]];then printf '\nsource "%s"\n' "$resource" >> "$destination/project.conf";fi
  sha256sum "$destination"/*.mcool > "$destination/input-sha256.txt"
fi
sha256sum --check "$destination/input-sha256.txt"
"$LOOP_ENV_PREFIX/bin/python" -I -B "$root/scripts/caller_probe.py" --mustache --chromosight
args=();if [[ "$resume" == true ]]; then args=(--resume);fi
bash "$root/scripts/microc2compare.sh" -c "$destination/project.conf" "${args[@]}"
run="$destination/results/demo/synthetic_demo/run001"
"$CORE_ENV_PREFIX/bin/python" -I -B "$root/scripts/validate_outputs.py" --config-json "$run/00_metadata/config.json"
bash "$root/scripts/microc2compare.sh" -c "$destination/project.conf" --resume
"$CORE_ENV_PREFIX/bin/python" -I -B - "$run" <<'PY'
import sys
from pathlib import Path
import numpy as np
import pandas as pd
root=Path(sys.argv[1])
for res in (2000,4000):
    counts=pd.read_csv(root/f'04_quantification/independent_loops/{res}/counts.tsv',sep='\t')
    assert len(counts)==1 and counts[['C1','C2','C3','T1','T2','T3']].notna().all().all()
    for sid in ('C1','C2','C3','T1','T2','T3'):
        scores=pd.read_csv(root/f'04_quantification/chromosight/{res}/{sid}.tsv',sep='\t')
        # Discovered-universe scoring may be empty; prespecified-only selection below must measure the fixture.
        assert set(scores.status)<= {'measured','ineligible_candidate','masked_window','flat_window','edge_window','undefined'}
    assert (root/f'05_aggregates/apa/{res}/C1.npz').is_file()
assert (root/'03_features/annotation/resources.json').is_file()
assert 'APA aggregates' in (root/'08_reports/index.html').read_text()
assert 'Compartment saddle plots' in (root/'08_reports/index.html').read_text()
partition=pd.read_csv(root/'05_aggregates/saddle/partition.tsv',sep='\t')
assert partition.partition.ge(0).sum()>=50, 'Saddle fixture lacks enough consistently oriented bins'
for unit in ('C1','C2','C3','T1','T2','T3','condition__control','condition__treated'):
    with np.load(root/f'05_aggregates/saddle/{unit}.npz') as data:
        assert np.isfinite(data['mean']).any(), f'No usable saddle values: {unit}'
    for suffix in ('png','svg'):
        assert (root/f'08_reports/figures/saddle.10000.{unit}.{suffix}').is_file()
    assert f'saddle.10000.{unit}.png' in (root/'08_reports/index.html').read_text()
print('Full fine-grid numerical products verified')
PY
# Independent-candidate selection tests scoring with no discovery requirement.
selected="$destination/selected.conf"
if [[ ! -e "$selected" ]]; then cp "$destination/project.conf" "$selected"; printf '\nRUN_ID="fixed_selected"\n' >> "$selected"; fi
selected_args=();if [[ -d "$destination/results/demo/synthetic_demo/fixed_selected" ]];then selected_args=(--resume);fi
bash "$root/scripts/microc2compare.sh" -c "$selected" --analysis chromosight,annotation,apa "${selected_args[@]}"
"$CORE_ENV_PREFIX/bin/python" -I -B - "$destination/results/demo/synthetic_demo/fixed_selected" <<'PY'
import sys
from pathlib import Path
import pandas as pd
root=Path(sys.argv[1])
for res in (2000,4000):
    for sid in ('C1','C2','C3','T1','T2','T3'):
        scores=pd.read_csv(root/f'04_quantification/chromosight/{res}/{sid}.tsv',sep='\t')
        assert len(scores)==1 and scores.status.eq('measured').all()
assert 'promoter_enhancer_like' in (root/'03_features/annotation/loops.tsv').read_text()
print('Selected fixed-candidate scoring, annotation and APA verified')
PY
# Saddle-only execution must also produce its report without loop discovery.
saddle="$destination/saddle.conf"
if [[ ! -e "$saddle" ]]; then cp "$destination/project.conf" "$saddle"; printf '\nRUN_ID="saddle_selected"\n' >> "$saddle"; fi
saddle_args=();if [[ -d "$destination/results/demo/synthetic_demo/saddle_selected" ]];then saddle_args=(--resume);fi
bash "$root/scripts/microc2compare.sh" -c "$saddle" --analysis saddle "${saddle_args[@]}"
"$CORE_ENV_PREFIX/bin/python" -I -B "$root/scripts/validate_outputs.py" --config-json "$destination/results/demo/synthetic_demo/saddle_selected/00_metadata/config.json"
"$CORE_ENV_PREFIX/bin/python" -I -B - "$run" "$destination/results/demo/synthetic_demo/saddle_selected" <<'PY'
import sys
from pathlib import Path
import numpy as np
full,selected=map(Path,sys.argv[1:])
for unit in ('C1','C2','C3','T1','T2','T3','condition__control','condition__treated'):
    with np.load(full/f'05_aggregates/saddle/{unit}.npz') as a, np.load(selected/f'05_aggregates/saddle/{unit}.npz') as b:
        np.testing.assert_allclose(a['mean'],b['mean'],rtol=1e-8,atol=1e-10,equal_nan=True)
    assert f'saddle.10000.{unit}.png' in (selected/'08_reports/index.html').read_text()
assert not (selected/'03_features/loops').exists()
print('Full/selected saddle arrays agree; all replicate and condition plots embedded in the final reports')
PY
sha256sum --check "$destination/input-sha256.txt"
printf 'Fine-loop acceptance passed. Report: %s/08_reports/index.html\n' "$run"
