#!/usr/bin/env bash
set -Eeuo pipefail
root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
while IFS= read -r -d '' file; do bash -n "$file"; done < <(find "$root/scripts" -name '*.sh' -print0)
bash "$root/scripts/microc2compare.sh" --help
"${PYTHON_BOOTSTRAP:-python3}" -m compileall -q "$root/src" "$root/scripts" "$root/tests"
