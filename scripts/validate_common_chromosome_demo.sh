#!/usr/bin/env bash
# Real Linux scientific acceptance; all writes are to a NEW private fixture.
set -Eeuo pipefail
umask 077
root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
[[ $# == 2 ]] || { printf '%s\n' 'Usage: bash scripts/validate_common_chromosome_demo.sh SOFTWARE.conf NEW_PRIVATE_DIRECTORY' >&2; exit 2; }
site=$1
destination=$2
[[ -r "$site" && ! -e "$destination" ]] || { printf '%s\n' 'Need a readable software profile and a new private destination' >&2; exit 2; }
source "$site"
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
"$CONDA_EXE" run --no-capture-output -p "$CORE_ENV_PREFIX" python "$root/scripts/generate_common_chromosome_demo.py" --out "$destination" --core-prefix "$CORE_ENV_PREFIX"
printf '\nsource "%s"\n' "$site" >> "$destination/project.conf"
sha256sum "$destination"/*.mcool > "$destination/input-sha256.txt"
bash "$root/scripts/microc2compare.sh" -c "$destination/project.conf" --dry-run
bash "$root/scripts/microc2compare.sh" -c "$destination/project.conf"
run="$destination/results/demo/common_chromosomes/run001"
"$PYTHON_BOOTSTRAP" "$root/scripts/validate_outputs.py" --config-json "$run/00_metadata/config.json"
"$CONDA_EXE" run --no-capture-output -p "$CORE_ENV_PREFIX" python - "$run" <<'PY'
import json,sys
from pathlib import Path
import numpy as np
import pandas as pd
root=Path(sys.argv[1])
for resolution in (10000,100000):
    base=root/f'00_metadata/masks/{resolution}'
    audit=json.loads(base.with_suffix('.scope.json').read_text())
    assert audit['analysis_chromosomes']==['chr1','chrX','chrY']
    assert audit['native_chromosome_sets_or_orders_differ']
    assert all(row['excluded_chromosomes']==(['scaffold'] if row['sample_id'].startswith('C') else []) for row in audit['samples'])
    assert audit['common_eligible_bins']==audit['analysis_grid_bins']
bins=pd.read_csv(root/'00_metadata/masks/10000.bins.tsv',sep='\t')
for sid in ('C1','C2','C3','T1','T2','T3'):
    table=pd.read_csv(root/f'03_features/insulation/{sid}.tsv',sep='\t')
    assert table[['chrom','start','end']].equals(bins)
    vectors=pd.read_csv(root/f'03_features/compartments/{sid}.tsv',sep='\t')
    assert set(vectors.chrom)=={'chr1','chrX','chrY'}
    assert set(pd.read_csv(root/f'01_qc/ps/{sid}.10000.tsv',sep='\t').region1)=={'chr1','chrX','chrY'}
pairwise=pd.read_csv(root/'01_qc/hicrep/by_view.tsv',sep='\t')
assert set(pairwise.region)=={'chr1','chrX','chrY'} and len(pairwise)==45
assert (root/'08_reports/index.html').is_file()
assert 'Chromosome scope and inherited weights' in (root/'08_reports/index.html').read_text()
assert not list(root.rglob('*.cool')) and not list(root.rglob('*.mcool'))
print('Full scientific common-chromosome/X/Y acceptance passed; no persistent derived matrices')
PY
# Selected module + resume exercise the same registration/coordinate mapping.
cp "$destination/project.conf" "$destination/tads-only.conf"
printf '\nRUN_ID="tads-only"\n' >> "$destination/tads-only.conf"
bash "$root/scripts/microc2compare.sh" -c "$destination/tads-only.conf" --analysis tads
bash "$root/scripts/microc2compare.sh" -c "$destination/tads-only.conf" --analysis tads --resume
"$PYTHON_BOOTSTRAP" "$root/scripts/validate_outputs.py" --config-json "$destination/results/demo/common_chromosomes/tads-only/00_metadata/config.json"
sha256sum --check "$destination/input-sha256.txt"
printf 'Common chromosome acceptance passed. Report: %s/08_reports/index.html\n' "$run"
