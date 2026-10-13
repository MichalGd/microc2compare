# Install and accept the 0.3.1 fine-loop candidate

Install alongside the active 0.2.x release. Existing analyses use their original versioned command. This upgrade reuses checked core, HiCRep and R prefixes; only a dedicated Mustache/Chromosight environment is new unless a compatible shared prefix already passes the numerical/API checks. The caller environment specification is unchanged from 0.3.0: reuse that prefix if it was already installed and passes the probes. The active shared alias changes only after acceptance. No user needs to modify shell startup files.

## Transfer and install

From **Windows PowerShell** transfer the generated archive, its checksum and installer from the workspace `.local` directory:

```powershell
$m2c_local = 'C:\work\microc2compare\.local'
scp -P PORT "$m2c_local\microc2compare-0.3.1-required-saddles-10oct26.zip" "$m2c_local\microc2compare-0.3.1-required-saddles-10oct26.zip.sha256" "$m2c_local\install-microc2compare-0.3.1.py" 'USER@SERVER:/opt/conda_envs/microc2compare/incoming/'
```

On **Linux**, in your normal Bash login:

```bash
base="/opt/conda_envs/microc2compare"
cd "$base/incoming"
sha256sum --check microc2compare-0.3.1-required-saddles-10oct26.zip.sha256
read -r checksum archive < microc2compare-0.3.1-required-saddles-10oct26.zip.sha256
/opt/miniconda/bin/python install-microc2compare-0.3.1.py \
  --archive "$archive" --sha256 "$checksum" --base "$base" \
  --release-id fine-loops-saddle-20261010
src="$base/releases/fine-loops-saddle-20261010/microc2compare-0.3.1"
```

The installer verifies the complete manifest before writing a new release. It refuses an existing release ID. It does not alter environments, upstream matrices or the shared alias.

## Reuse core/R; create or verify only fine-loop software

Prepare a new maintainer-owned software profile by copying the accepted profile, and append the new loop prefix:

```bash
umask 022
site="$base/site/shared-fine-loops-saddle-20261010.conf"
test ! -e "$site" && {
cp "$base/site/shared-r-ql-fix-20261008.conf" "$site" &&
printf '\nLOOP_ENV_PREFIX="%s"\nLOOP_LOCK=""\n' \
  "$base/envs/fine-loops-0.3.0-20261010" >> "$site" &&
chmod 644 "$site"
}
```

Run this creation block once; if the profile exists, inspect it before editing rather than repeatedly appending lines. To reuse an already installed Mustache/Chromosight prefix, set its path instead and run `caller_probe.py` from its Python. Require Mustache 1.3.3, Chromosight 1.6.3 and the checked API; a package name alone is insufficient.

For a new prefix, keep installation caches in your own workspace and software readable to all users:

```bash
work="$HOME/microc2compare-work"
mkdir -p "$work/cache/conda-pkgs" "$work/cache/pip" "$work/cache/matplotlib" "$work/tmp" "$work/audit"
# If this exact prefix already has conda-meta, reuse it and run the probe below.
test -d "$base/envs/fine-loops-0.3.0-20261010/conda-meta" || {
test ! -e "$base/envs/fine-loops-0.3.0-20261010" &&
CONDA_REGISTER_ENVS=false CONDA_ALWAYS_COPY=true \
CONDA_PKGS_DIRS="$work/cache/conda-pkgs" PIP_CACHE_DIR="$work/cache/pip" TMPDIR="$work/tmp" \
/opt/miniconda/bin/conda env create \
  -p "$base/envs/fine-loops-0.3.0-20261010" -f "$src/envs/fine-loops.yml"
}
MPLCONFIGDIR="$work/cache/matplotlib" "$base/envs/fine-loops-0.3.0-20261010/bin/python" -I -B \
  "$src/scripts/caller_probe.py" --mustache --chromosight
```

No package update is applied to the reused core/R/HiCRep environments. The YAML is a versioned dependency specification, not an already tested exact Linux lock. After acceptance record `conda list --explicit`, `conda list --json` and `pip freeze`, including pip-installed caller packages. Exact reconstruction of a mixed environment needs reviewed hashed pip requirements as well as the Conda explicit list.

## Candidate acceptance before activation

Publish a versioned command, binding the profile internally:

```bash
/opt/miniconda/bin/python "$src/scripts/publish_shared_command.py" \
  --source "$src" --base "$base" --site-profile "$site"
candidate="$base/bin/microc2compare-0.3.1-fine-loops-saddle-20261010"
"$candidate" --list-analyses
```

Run the dedicated synthetic acceptance in a **new private** directory. It exercises actual Linux caller imports/discovery, complete 2/4-kb measurement, fixed Chromosight scoring, annotation, APA/ADA, required compartment/saddle analysis and embedded replicate/condition saddle plots, descriptive effects, selected modules, resume and unchanged input hashes. It is software validation, not biological calibration.

```bash
umask 077
check="$HOME/microc2compare-work/validation/fine-loops-031-001"
nohup bash "$src/scripts/validate_fine_loop_demo.sh" \
  --site "$site" --out "$check" \
  > "$HOME/microc2compare-work/audit/fine-loops-031-001.log" 2>&1 < /dev/null &
tail -f "$HOME/microc2compare-work/audit/fine-loops-031-001.log"
```

Ctrl+C exits `tail` without stopping the background job. Reconnect and inspect the log to confirm the final acceptance marker. To continue an interrupted acceptance use the same command with `--resume`; do not start a second process on the same run.

Check non-owner command/prefix access before activation:

```bash
sudo -u nobody env -i PATH=/usr/bin:/bin HOME=/nonexistent \
  "$candidate" --list-analyses
sudo -u nobody env -i PATH=/usr/bin:/bin HOME=/nonexistent \
  "$base/envs/fine-loops-0.3.0-20261010/bin/python" -I -B \
  "$src/scripts/caller_probe.py" --mustache --chromosight
```

Shared source/profile/bin/prefix directories should be readable/executable and not writable by ordinary researchers. Basic import checks do not replace a full run by a second research account with its own input/output/scratch directories. The launcher internally sets a private umask, clears Python/R contamination and places caches under the configured scratch tree; ordinary users run `microc2compare -c project.conf`.

After numerical/server acceptance, activate explicitly:

```bash
/opt/miniconda/bin/python "$src/scripts/publish_shared_command.py" \
  --source "$src" --base "$base" --site-profile "$site" --activate
hash -r
readlink -f "$base/bin/microc2compare"
```

Keep the older versioned launcher for historical runs and rollback. A real-data pilot and scope/resolution-specific statistical calibration remain required before formal biological inference. The local verification record lists what was actually tested; neither installation nor a synthetic demonstration establishes calibration.
