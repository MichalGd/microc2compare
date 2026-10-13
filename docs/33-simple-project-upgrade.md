# Install and accept the 0.3.2 convenience candidate

This is a **maintainer procedure**, performed once for a shared installation. Researchers subsequently use [simple project commands](32-simple-project-commands.md). The candidate changes project interfaces/documentation and preserves the fine-loop scientific methods, mandatory saddle reports and statistical gates from 0.3.1. Core, HiCRep, R and optional caller dependency definitions are unchanged. Reuse accepted shared environments; do not update other workflows' environments.

## Transfer and install a separate immutable release

In **Windows PowerShell**:

```powershell
$local = 'C:\work\microc2compare\.local'
scp -P PORT "$local\microc2compare-0.3.2-simple-projects-10oct26.zip" "$local\microc2compare-0.3.2-simple-projects-10oct26.zip.sha256" "$local\install-microc2compare-0.3.2.py" 'USER@SERVER:/opt/conda_envs/microc2compare/incoming/'
```

In your normal **Linux Bash** login:

```bash
base="/opt/conda_envs/microc2compare"
cd "$base/incoming"
sha256sum --check microc2compare-0.3.2-simple-projects-10oct26.zip.sha256
read -r checksum archive < microc2compare-0.3.2-simple-projects-10oct26.zip.sha256
archive="${archive%$'\r'}"
/opt/miniconda/bin/python install-microc2compare-0.3.2.py \
  --archive "$archive" --sha256 "$checksum" --base "$base" \
  --release-id simple-projects-20261010 &&
src="$base/releases/simple-projects-20261010/microc2compare-0.3.2"
```

The `archive` line removes a trailing Windows carriage return if the checksum file uses CRLF line endings; it also works with Linux LF files. If the installer reports a filename ending in `\r`, repeat this corrected block. That filename error occurs before the installer creates the release directory, so the same release ID can be reused.

The installer verifies the manifest and refuses an existing release ID. It preserves the active shared command, old releases, environments and project data.

## Bind checked software and test

Use a software profile including the separate Mustache/Chromosight prefix for the full fine-loop strategy. **The release installer does not create site profiles or environments.** The path below is the proposed new profile, and it must exist before publication. Your previously accepted core, HiCRep and R profile can be reused to create it once:

```bash
site="$base/site/shared-fine-loops-saddle-20261010.conf"
if [[ ! -e "$site" ]]; then
  (
    set -e
    umask 022
    test -r "$base/site/shared-r-ql-fix-20261008.conf"
    cp "$base/site/shared-r-ql-fix-20261008.conf" "$site"
    printf '\nLOOP_ENV_PREFIX="%s"\nLOOP_LOCK=""\n' \
      "$base/envs/fine-loops-0.3.0-20261010" >> "$site"
    chmod 644 "$site"
  )
fi
```

An existing profile is preserved. If a compatible checked caller prefix has already been registered elsewhere, use that accepted profile instead. The 0.3.0 caller definition remains compatible with this 0.3.2 release. These are one-time maintainer steps, not normal per-project setup. Creating a profile or displaying command help does not establish that its environments work.

If `$base/envs/fine-loops-0.3.0-20261010/bin/python` is missing, create the dedicated caller environment using the definition in the **current 0.3.2 release**. The following command survives an SSH disconnect, keeps installation caches/logs in your home workspace, and sets a readable software umask in the background process. It reuses an existing Conda prefix for the probe and refuses to overwrite a partial directory. Core, HiCRep and R environments are reused as they are.

```bash
base="/opt/conda_envs/microc2compare"
src="$base/releases/simple-projects-20261010/microc2compare-0.3.2"
loops="$base/envs/fine-loops-0.3.0-20261010"
work="$HOME/microc2compare-work"
umask 077
mkdir -p "$work"/{audit,tmp,cache/conda-pkgs,cache/pip,cache/matplotlib}
log="$work/audit/fine-loops-install-$(date -u +%Y%m%d-%H%M%S).log"

nohup bash -c '
set -euo pipefail
umask 022
unset PYTHONHOME PYTHONPATH
prefix=$1
src=$2
work=$3
export CONDA_REGISTER_ENVS=false CONDA_ALWAYS_COPY=true CONDA_ALWAYS_YES=true
export CONDA_PKGS_DIRS="$work/cache/conda-pkgs"
export PIP_CACHE_DIR="$work/cache/pip" TMPDIR="$work/tmp"
export MPLCONFIGDIR="$work/cache/matplotlib" PYTHONNOUSERSITE=1
test -r "$src/envs/fine-loops.yml"
if [[ ! -d "$prefix/conda-meta" ]]; then
  test ! -e "$prefix" || {
    printf "ERROR: existing incomplete prefix; inspect it before retrying: %s\n" "$prefix" >&2
    exit 2
  }
  /opt/miniconda/bin/conda env create \
    -p "$prefix" -f "$src/envs/fine-loops.yml"
fi
"$prefix/bin/python" -I -B "$src/scripts/caller_probe.py" --mustache --chromosight
printf "Fine-loop environment setup and probe passed\n"
' _ "$loops" "$src" "$work" > "$log" 2>&1 < /dev/null &

printf 'PID: %s\nLog: %s\n' "$!" "$log"
tail -f "$log"
```

Ctrl+C exits `tail` only. Wait for `Fine-loop environment setup and probe passed` before proceeding. On reconnection, restore `base` and `src` from the first two lines above and inspect the printed log path with `tail -n 80 PATH`. A successful package/API probe is followed by both acceptance fixtures below; it is not full workflow acceptance. The [fine-loop guide](31-fine-loop-upgrade.md#reuse-corer-create-or-verify-only-fine-loop-software) also describes recording the accepted dependency inventory.

The publisher error `Shared software profile must be a file below the shared site directory` means the resolved path is absent, is not a regular file, or falls outside `site`. With the illustrated path, check `ls -l "$site"` and complete the creation step above. Publication stops before writing the candidate launcher, so subsequent `No such file or directory` errors for that launcher are expected. The active command remains unchanged.

Once the profile exists, publish the versioned command without activating it. Keep the commands connected with `&&` so help checks run only after successful publication:

```bash
test -r "$site" &&
/opt/miniconda/bin/python "$src/scripts/publish_shared_command.py" \
  --source "$src" --base "$base" --site-profile "$site" &&
candidate="$base/bin/microc2compare-0.3.2-simple-projects-20261010" &&
"$candidate" --help &&
"$candidate" --list-analyses
```

Before shared activation run the convenience acceptance in a new private directory. It creates a compact configuration, generates the known synthetic QC fixture using the installed core, audits its resolutions, executes QC/HiCRep through the bound command, validates completion and confirms unchanged maps. It needs the accepted core and HiCRep prefixes and the small fixture's 4-CPU/16-GB budget.

```bash
work="$HOME/microc2compare-work"
mkdir -p "$work/audit"
check="$work/validation/simple-projects-032-001"
nohup bash "$src/scripts/validate_project_tools_demo.sh" \
  --command "$candidate" --out "$check" \
  > "$work/audit/simple-projects-032-001.log" 2>&1 < /dev/null &
tail -f "$work/audit/simple-projects-032-001.log"
```

This is software validation, not statistical calibration. Also run the full fine-loop acceptance from [the fine-loop guide](31-fine-loop-upgrade.md#candidate-acceptance-before-activation), using this `src`, the accepted `site` and a new private destination. That gate executes actual callers, compartments, saddles, reports, selected modules and resume; local Windows checks cannot replace it.

## Ordinary-account access and activation

After both owner-run acceptance fixtures pass, verify shared-file permissions and actual workflow execution as a non-owner. The command below uses `nobody` and a private temporary directory containing only synthetic data. It does not require access to your home directory or create a shared research-results folder. The temporary audit directory is retained and printed at completion.

```bash
sudo -u nobody env -i PATH=/usr/bin:/bin HOME=/nonexistent \
  bash --noprofile --norc -c '
set -euo pipefail
umask 077
cd /
base=$1
src=$2
site=$3
candidate=$4
source "$site"
for path in "$base" "$base/site" "$base/bin" "$src" "$site" \
  "$CORE_ENV_PREFIX" "$HICREP_ENV_PREFIX" "$DIFF_R_ENV_PREFIX" "$LOOP_ENV_PREFIX"; do
  test -r "$path"
  test ! -w "$path"
done
check=$(mktemp -d /tmp/microc2compare-032-access.XXXXXX)
mkdir -p "$check"/{tmp,cache,matplotlib}
export HOME="$check" TMPDIR="$check/tmp"
export XDG_CACHE_HOME="$check/cache" MPLCONFIGDIR="$check/matplotlib"
"$LOOP_ENV_PREFIX/bin/python" -I -B "$src/scripts/caller_probe.py" --mustache --chromosight
bash "$src/scripts/validate_project_tools_demo.sh" \
  --command "$candidate" --out "$check/demo"
printf "Non-owner workflow execution passed. Audit: %s\n" "$check"
' _ "$base" "$src" "$site" "$candidate"
```

The same convenience fixture can instead be run under another researcher's ordinary account in **their own writable project directory**. Shared code, profiles and prefixes must be readable/executable but not writable by ordinary users. A help check alone is not a full user execution check. No `.bashrc` change, Conda initialization or shared writable results area is needed.

After acceptance:

```bash
/opt/miniconda/bin/python "$src/scripts/publish_shared_command.py" \
  --source "$src" --base "$base" --site-profile "$site" --activate
hash -r
readlink -f "$base/bin/microc2compare"
microc2compare --help
```

Existing runs should continue with their original versioned launcher/configuration. For a new fine-grid project use a new run identity and inspect genome/phasing, coverage, reproducibility and resources before full real-data analysis. Formal biological inference still requires applicable calibration and design gates. This local delivery does not itself install, activate or execute anything on the server.
