#!/usr/bin/env bash
# One-time maintainer/ordinary-account Linux acceptance; no shared writes.
set -Eeuo pipefail
command= out=
while (( $# )); do
  case "$1" in
    --command|--out)
      (( $# >= 2 )) || { printf 'Missing value for %s\n' "$1" >&2; exit 2; }
      case "$1" in --command) command=$2;; --out) out=$2;; esac
      shift 2;;
    *) printf 'Unknown option: %s\n' "$1" >&2; exit 2;;
  esac
done
[[ -n "$command" && -x "$command" && -n "$out" && ! -e "$out" ]] || { printf 'Use --command CHECKED_COMMAND --out NEW_PRIVATE_DIRECTORY\n' >&2; exit 2; }
umask 077
"$command" --init "$out/starter/project.conf"
# Resolve the accepted software through the candidate itself, without asking
# users to provide fixed interpreter or source paths.
"$command" -c "$out/starter/project.conf" --demo --out "$out/data"
conf="$out/data/project.conf"
find "$out/data" -maxdepth 1 -name '*.mcool' -print0 | sort -z | xargs -0 sha256sum > "$out/inputs.sha256"
"$command" -c "$conf" --audit-resolutions --out "$out/audit"
"$command" -c "$conf" --analysis qc,hicrep
"$command" -c "$conf" --validate
sha256sum --check "$out/inputs.sha256"
test -s "$out/data/results/demo/synthetic_demo/run001/08_reports/index.html"
printf 'Project command acceptance passed; software selection and private outputs verified\n'
