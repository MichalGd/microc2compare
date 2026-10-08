# Install once for all users: shared clone, separate workspaces

For a new installation from repository source, use [the GitHub installation guide](26-install-from-github.md). ZIP transfer sections here document earlier private deliveries; those artifacts are not included in Git.

**Recommended deployment:** one maintainer-owned application under `/opt/conda_envs/microc2compare`, one versioned Conda clone of the existing Linux `microc2tracks` core, and a checked site profile/command for every user. Reuse an existing centrally accessible Conda manager outside personal homes. Make any necessary package additions only to the new clone. Create a separate HiCRep environment only if needed; defer differential R until formal analysis is prepared. Keep the source environment, your home installation, shell startup files and existing workflows unchanged.

Keep each user's installation audits, package caches, temporary files and analysis projects in their own home directory, by default `$HOME/microc2compare-work`. This new workspace is separate from existing projects and from the shared software tree. Users may choose another personal project directory, but `/dysk3` is not required. Ordinary users run the published software without cloning/installing it themselves. This deliberate one-time environment copy costs disk space but separates future `microc2compare` changes from upstream work. Raw counts/stored weights and matching references are still read in place; cloning software does not request copies or rebalancing of upstream matrices.

This guide prepares a shared **exploratory** installation. Complete [publication and ordinary-account acceptance](16-shared-server-installation.md) before announcing access to all users. Linux/R/calibration/public-data/server validation remains pending until actually executed and recorded.

Run the steps in order as the designated maintainer. Commands marked **Windows PowerShell** run on your workstation; all other commands run in a **Bash shell on the Linux server**. Replace `USERNAME`, `SERVER` and `/ACTUAL/...` with real values. Stop at a failed check and read its log before proceeding. Use a dedicated terminal for installation; the session-only exports below do not alter your other terminals. Keep the Linux terminal open so the `m2c_*` variables remain available; after reconnecting, re-enter the path/selection assignments from steps 1, 2, 3 and 6. This guide's `$m2c_base` contains shared software; `$m2c_work` contains this account's work.

## 1. Log in and select shared software and private work directories

**Windows PowerShell:**

```powershell
ssh USERNAME@SERVER
```

**Linux server:** use an administrator-approved, maintainer-owned shared software directory, such as `/opt/conda_envs/microc2compare`. It must be reachable by every intended account through all parent directories. If it has not been provisioned, first use [the administrator preparation step](16-shared-server-installation.md#1-provision-the-shared-software-root-once). Do not install the shared software inside a personal home or open a home directory to other users. Use another shared application filesystem only if site policy requires it. Use `$HOME/microc2compare-work` for your own work; other accounts use their own `$HOME`. No workspace on `/dysk3`, administrator-provisioned project directory or home permission change is required.

Run this **single command first**, press Enter, and wait for the new Bash prompt before pasting anything else:

```bash
bash --noprofile --norc
```

Then run the next block. `$HOME` expands to the current account's home directory; it is not a placeholder requiring replacement. The directory checks give an explanation if the shared root has not yet been provisioned. A new Bash prompt alone does not mean installation succeeded; look for the printed `Shared software` and `Private work` paths and filesystem output.

```bash
set -Eeuo pipefail
unset HISTFILE BASH_ENV ENV PYTHONHOME PYTHONPATH R_HOME R_LIBS R_LIBS_USER
m2c_base="/opt/conda_envs/microc2compare"
m2c_work="$HOME/microc2compare-work"
if [[ ! -d "$m2c_base" || ! -w "$m2c_base" ]]; then
  printf 'Shared root missing or not writable: %s\nUse shared-guide step 1 to provision it first.\n' "$m2c_base" >&2
  exit 2
fi
umask 022
mkdir -p "$m2c_base/incoming" "$m2c_base/releases" "$m2c_base/envs" \
  "$m2c_base/site" "$m2c_base/bin"
umask 077
mkdir -p "$m2c_work/projects" "$m2c_work/audit" "$m2c_work/validation" \
  "$m2c_work/tmp" "$m2c_work/cache/conda-pkgs" \
  "$m2c_work/runtime-config" "$m2c_work/runtime-data" \
  "$m2c_work/cache/matplotlib" "$m2c_work/cache/pip"
export TMPDIR="$m2c_work/tmp"
export CONDA_PKGS_DIRS="$m2c_work/cache/conda-pkgs"
export CONDA_ENVS_PATH="$m2c_base/envs"
export CONDA_REGISTER_ENVS=false
export CONDA_ALWAYS_COPY=true
export XDG_CACHE_HOME="$m2c_work/cache"
export XDG_CONFIG_HOME="$m2c_work/runtime-config"
export XDG_DATA_HOME="$m2c_work/runtime-data"
export MPLCONFIGDIR="$m2c_work/cache/matplotlib"
export PIP_CACHE_DIR="$m2c_work/cache/pip"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
printf 'Shared software: %s\nPrivate work: %s\n' "$m2c_base" "$m2c_work"
df -hT "$m2c_base"
df -hT "$m2c_work"
df -i "$m2c_base"
if command -v quota >/dev/null 2>&1; then
  quota -s || printf '%s\n' 'Quota information unavailable; confirm capacity with the site.'
fi
bash --version
uname -m
```

Expected: maintainer-writable software root, private writable work directories, sufficient disk/inodes/quota, Linux `x86_64`, Bash 4.4+. `quota` may be unavailable or report no quota; record that and confirm filesystem capacity by the server's normal procedure. Administrator involvement is needed only for site directory/global command provisioning when the maintainer lacks that authority; ordinary analysis never needs sudo. Use the site's approved compute allocation for tests if the server separates login and compute nodes. The private package cache avoids writes to an upstream cache and needs space for cloning/downloads. `CONDA_ALWAYS_COPY` keeps newly provisioned package files separate from caches; `CONDA_REGISTER_ENVS=false` avoids adding the new prefixes to your home environment registry. Unsetting `HISTFILE` prevents this dedicated shell from saving command history in home; clearing `BASH_ENV`/`ENV` avoids inherited shell startup hooks in child commands. All settings are confined to this terminal and its children. Do not run `conda init`, `conda config --set`, `conda update`, or change home/upstream permissions. Conda registration control is version-dependent and is checked in step 2; do not bypass that check.

## 2. Select the existing shared Conda manager

Use the manager already used by the shared RNA/ATAC installation, after confirming its actual location and access. The RNA repository documents `/opt/miniconda`; this is a candidate, not a verified server fact. Use the executable directly, without activation or shell initialization:

```bash
m2c_conda="/opt/miniconda/bin/conda"
test -x "$m2c_conda"
"$m2c_conda" --version
m2c_conda_base=$("$m2c_conda" info --base)
printf 'Conda executable: %s\nConda base: %s\n' "$m2c_conda" "$m2c_conda_base"
case "$m2c_conda_base" in
  /home/*|/root|/root/*) printf '%s\n' 'Select a centrally accessible manager outside personal homes.' >&2; exit 2 ;;
esac
"$m2c_conda" config --show register_envs
```

Expected: a working centrally accessible manager and `register_envs: false`. Confirm the executable's resolved path and all parents are outside personal homes too. If `/opt/miniconda/bin/conda` does not exist, obtain the existing manager's actual path from the site installation record; do not install a duplicate by default. If the only manager is inside your home, the administrator must provide an independent central manager. Leave the home installation in place; do not move it or broaden its access.

The microc2tracks server guide documents core prefix `/opt/conda/envs/microc2tracks` but manager `/home/USERNAME/miniconda3`. Only the core is the cloning candidate; the home-based manager is not selected for this deployment. These are local repository records, not live server checks.

Record known shared prefixes from the administrator/installation records and inspect them by absolute path. Avoid `conda env list`/`conda info --envs`: some versions clean `~/.conda/environments.txt` while listing. No installation step needs those commands. Conda's [multi-user guidance](https://docs.conda.io/projects/conda/en/stable/user-guide/configuration/admin-multi-user-install.html) describes centrally managed installations; registration behavior is implemented in [Conda's environment manager](https://github.com/conda/conda/blob/main/conda/core/envs_manager.py). Verify this control on the actual installed version before creating a prefix.

Record selected home settings before installation, for comparison in step 9:

```bash
for m2c_home_file in "$HOME/.bashrc" "$HOME/.bash_profile" "$HOME/.profile" \
  "$HOME/.condarc" "$HOME/.conda/environments.txt"; do
  if [[ -f "$m2c_home_file" ]]; then
    sha256sum "$m2c_home_file"
  else
    printf 'ABSENT  %s\n' "$m2c_home_file"
  fi
done > "$m2c_work/audit/home-settings-before.txt"
```

This reads only the named settings; it does not modify them or prove every possible third-party hook/cache is free of home writes. Review activation hooks and run the acceptance checks below. Keep audits private because paths can identify accounts.

## 3. Select an existing Python 3.11+ for setup commands

The bootstrap Python runs small standard-library configuration/metadata commands. It does not have to be the scientific environment. Try the existing Conda base Python:

```bash
m2c_conda_base=$("$m2c_conda" info --base)
m2c_bootstrap="$m2c_conda_base/bin/python"
"$m2c_bootstrap" -c 'import sys; print(sys.executable); print(sys.version); assert sys.version_info >= (3, 11)'
```

If too old, try an existing system Python or another central prefix, then repeat the version check. Its resolved executable and dependencies must be accessible outside personal homes:

```bash
m2c_bootstrap="$(type -P python3)"
"$m2c_bootstrap" -c 'import sys; print(sys.executable); print(sys.version); assert sys.version_info >= (3, 11)'
```

Only if **no existing Python 3.11+ accessible to all users is available**, use the existing Conda executable to create a small, new maintainer-owned shared bootstrap prefix. Apply the shared-software umask to this installation command:

```bash
m2c_bootstrap_prefix="$m2c_base/envs/bootstrap-0.1.0"
test ! -e "$m2c_bootstrap_prefix" && \
  (umask 022; "$m2c_conda" create -y -p "$m2c_bootstrap_prefix" -c conda-forge python=3.11)
m2c_bootstrap="$m2c_bootstrap_prefix/bin/python"
"$m2c_bootstrap" -c 'import sys; print(sys.version); assert sys.version_info >= (3, 11)'
```

If that prefix already exists, inspect it or choose a new name; the creation guard must not be bypassed. Do not upgrade Conda base or an upstream Python to satisfy this requirement.

## 4. Transfer and verify the reviewed source archive

The current locally prepared archive is `microc2compare-0.1.0-design-json-fix-8oct26.zip`. It includes the compartment export, masked bin-selector and design JSON corrections, and supersedes the earlier candidates. There is no published microc2compare repository assumed by these commands. If you already installed and checked the core and HiCRep prefixes, retain them; install the new source in a separate release directory and resume the existing fixture rather than cloning or installing the environments again. Full-fixture and ordinary-account acceptance still apply before command publication.

Open a **second Windows PowerShell terminal**; leave the Linux session open:

```powershell
scp "D:\path\to\microc2compare\.local\microc2compare-0.1.0-design-json-fix-8oct26.zip" "USERNAME@SERVER:/opt/conda_envs/microc2compare/incoming/"
```

This destination assumes the shared `/opt/conda_envs/microc2compare` from step 1. For a different base, replace the remote destination with its actual absolute `.../incoming/` path. Log in as the maintainer account that owns it.

**Linux server:**

```bash
m2c_archive="$m2c_base/incoming/microc2compare-0.1.0-design-json-fix-8oct26.zip"
test -f "$m2c_archive"
"$m2c_bootstrap" -m zipfile -t "$m2c_archive"
m2c_release="$m2c_base/releases/shared-design-json-fix-20261008"
test ! -e "$m2c_release" && \
  (umask 022; mkdir "$m2c_release" && \
    "$m2c_bootstrap" -m zipfile -e "$m2c_archive" "$m2c_release")
m2c_src="$m2c_release/microc2compare-0.1.0"
cd "$m2c_src"
```

Expected: archive integrity succeeds and the source folder contains `README.md`, `scripts/`, `config/` and `RELEASE_MANIFEST.json`. If the release directory exists, choose a new release-directory name before extraction. Keep older releases in place for rollback.

Verify the recorded source checksums:

```bash
"$m2c_bootstrap" - "$m2c_src" <<'PY'
import hashlib, json, sys
from pathlib import Path
root = Path(sys.argv[1])
manifest = json.loads((root / 'RELEASE_MANIFEST.json').read_text())['files']
for name, expected in manifest.items():
    assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, name
print(f'Verified {len(manifest)} source/documentation files')
PY
bash scripts/microc2compare.sh --help
"$m2c_bootstrap" scripts/prepare_microc2compare_inputs.py --help
```

Expected: every source checksum matches, and both help commands succeed. No package installation into an existing environment is needed merely to run from this checkout.

## 5. Record the real server allocation and storage

```bash
PATH="$(dirname "$m2c_conda"):$PATH" \
  bash scripts/server_inventory.sh > "$m2c_work/audit/server-inventory.txt" 2>&1
cat "$m2c_work/audit/server-inventory.txt"
```

Review CPU affinity/quota, memory limits/current availability, filesystem space/quota, CPU/NUMA topology and existing workload. The reported 500 GB RAM, four Xeons, 70 physical cores and 140 logical threads are a hardware description, not your automatic allocation. Missing optional inventory utilities are recorded; the script does not install them or change system settings.

| Step | Initial budget | Purpose |
|---|---|---|
| Synthetic smoke/full-feature demonstration | 4 CPUs, 16,000 MB; one sample/matrix worker | Verify installation with tiny data |
| First representative real-data pilot | Up to 24 CPUs, 160,000 MB; two sample workers, one I/O-heavy worker | Measure memory, runtime and storage behavior |
| Later production candidate | Up to 60 CPUs, 380,000 MB | Consider only after inventory/pilot and available allocation permit it |

Memory units are decimal MB. These are admission budgets, not hard operating-system memory limits. Do not allocate 140 CPU workers because hyperthreading reports 140 logical threads. Use [server-resource guidance](07-server-resources.md) and the site's scheduler/allocation policy.

For a smaller permitted allocation, lower both global budgets and enabled task reservations before running the demo. Stop if the current runtime reports unsupported cgroup v1 or insufficient allocation; do not change the server's kernel/cgroup setup merely to pass preflight.

## 6. Clone the microc2tracks core and check the clone

Use the actual installed **Linux core**, not the Windows VCS folder or the optional downstream environment. The local core definition declares Python 3.11, Cooler, cooltools, bioframe, NumPy, SciPy and pandas; most versions are unpinned. It does not explicitly declare `jsonschema` or HiCRep. Inventory the real prefix before deciding which additions are needed.

```bash
m2c_upstream_core="/opt/conda/envs/microc2tracks"
m2c_core="$m2c_base/envs/core-0.1.0-clone-20261007"
test -d "$m2c_upstream_core/conda-meta"
test ! -e "$m2c_core"
"$m2c_conda" list -p "$m2c_upstream_core" --json \
  > "$m2c_work/audit/upstream-core-packages-before.json"
find "$m2c_upstream_core/conda-meta" -maxdepth 1 -type f -print0 | \
  sort -z | xargs -0 -r sha256sum > "$m2c_work/audit/upstream-conda-meta-before.txt"
du -sh "$m2c_upstream_core"
df -h "$m2c_base/envs" "$m2c_work/cache/conda-pkgs"
"$m2c_conda" run --no-capture-output -p "$m2c_upstream_core" \
  python -c 'import sys; print(sys.version); assert sys.version_info >= (3, 11)'
(umask 022; "$m2c_conda" create -y --no-default-packages \
  --prefix "$m2c_core" --clone "$m2c_upstream_core" --copy) \
  2>&1 | tee "$m2c_work/audit/core-clone.log"
```

Confirm capacity for both separate environment files and the private package cache before cloning. Downloads may still be needed; cloning is not guaranteed offline. Do not clone while the upstream maintainer is updating its environment. If the source is missing, inaccessible or incompatible at the Python check, skip cloning, retain the audit, create the software-check file below with the proposed new `m2c_core`, and use step 7. With fail-fast enabled, a failed command may close the installation subshell; restart it and restore the selections/session exports before following the fallback. An existing target is a stop point: inspect it or choose a new versioned name.

Conda [supports managed cloning and `--copy`](https://docs.conda.io/projects/conda/en/stable/commands/create.html). Ordinary filesystem copying/moving is not the installation method. `--copy` avoids package linking during creation but does not remove arbitrary external symlinks, editable pip installs or custom activation-hook dependencies. Inspect those dependencies and test the clone without a home-based manager; use the clean fallback in step 7 if it depends on private home software. Upstream FASTQ/matrices/reference directories are outside the environment and are not requested as part of this software clone.

Create a **software-only configuration** outside the release. It has no real sample tables. Differential analysis is disabled while core feature switches stay on to check the required scientific APIs:

```bash
m2c_setup="$m2c_work/projects/software-check"
mkdir -p "$m2c_setup"
m2c_setup_conf="$m2c_setup/setup.conf"
test ! -e "$m2c_setup_conf"
cat > "$m2c_setup_conf" <<EOF
source "$m2c_src/config/resources_demo.conf"
ENVIRONMENT_MODE="reuse"
CONDA_EXE="$m2c_conda"
PYTHON_BOOTSTRAP="$m2c_bootstrap"
CORE_ENV_PREFIX="$m2c_core"
HICREP_ENV_PREFIX="$m2c_core"
RUN_HICREP="false"
RUN_DIFFERENTIAL="false"
RUN_COMPARTMENTS="true"
RUN_INSULATION="true"
RUN_LOOPS="true"
RUN_POOL="false"
RUN_DEPTH_SENSITIVITY="false"
BALANCE_NEW_DERIVED_MAPS="false"
OUTDIR="$m2c_setup/results"
SCRATCH_ROOT="$m2c_setup/scratch"
EOF
```

If cloning succeeded, check it. A failed probe is expected only when diagnosing a candidate; the conditional distinguishes success from failure without terminating the installation shell:

```bash
if bash scripts/bootstrap_environments.sh -c "$m2c_setup_conf" --check-only \
  2>&1 | tee "$m2c_work/audit/core-check.log"; then
  printf '%s\n' 'Clone passed the core probe; skip step 7.'
else
  printf '%s\n' 'Clone failed: review the log, then use the correction or fallback below.'
fi
```

Expected: exit 0, numerical raw/weighted query checks pass and core feature APIs are available. The probe creates a tiny disposable map, never changes real upstream matrices, and is not a substitute for step 10's full Linux fixture.

If the **only** reported deficiency is missing `jsonschema`, review a target-only dry solve first:

```bash
"$m2c_conda" install -p "$m2c_core" --override-channels \
  -c conda-forge -c bioconda --freeze-installed --dry-run jsonschema=4.23.0 \
  > "$m2c_work/audit/clone-jsonschema-dry-run.txt"
cat "$m2c_work/audit/clone-jsonschema-dry-run.txt"
```

Only when that plan leaves the existing scientific packages unchanged, run the matching installation into the clone and recheck:

```bash
(umask 022; "$m2c_conda" install -y -p "$m2c_core" --override-channels \
  -c conda-forge -c bioconda --freeze-installed jsonschema=4.23.0) \
  2>&1 | tee "$m2c_work/audit/clone-jsonschema-install.log"
bash scripts/bootstrap_environments.sh -c "$m2c_setup_conf" --check-only \
  2>&1 | tee "$m2c_work/audit/core-recheck.log"
```

Do not execute this correction for an absent/failed clone, other missing APIs, NumPy incompatibility or a solver plan changing the scientific stack. Use the fresh pinned fallback instead. No command installs packages into `m2c_upstream_core` or upgrades Conda base. Once accepted/published, the clone is held unchanged; later additions require a new environment/profile version.

## 7. Clean fallback only if cloning/checking failed

**Skip when the clone passed.** Otherwise select a fresh shared destination; keep the failed clone/logs for diagnosis and do not overwrite or automatically delete them:

```bash
m2c_core="$m2c_base/envs/core-0.1.0-clean-20261007"
test ! -e "$m2c_core"
printf 'CORE_ENV_PREFIX="%s"\n' "$m2c_core" >> "$m2c_setup_conf"
(umask 022; bash scripts/bootstrap_environments.sh -c "$m2c_setup_conf" --create-missing) \
  2>&1 | tee "$m2c_work/audit/core-create.log"
bash scripts/bootstrap_environments.sh -c "$m2c_setup_conf" --check-only \
  2>&1 | tee "$m2c_work/audit/core-recheck.log"
"$m2c_conda" run --no-capture-output -p "$m2c_core" \
  python -c 'import sys; print(sys.version); assert sys.version_info >= (3, 11)'
```

This uses the existing central manager and checked-in `environment.yml` to create only a new `microc2compare` core. It does not alter existing Conda prefixes. Existing non-Conda destinations are refused. Preserve solve/probe logs and choose a fresh suffix if a previous attempt failed. Principal package pins are candidates, not a fully validated Linux binary lock. Expected: creation and both checks pass; unresolved failures stop installation.

## 8. Reuse HiCRep if available; otherwise isolate it

Try the successful core prefix first, or another **existing** prefix known to contain HiCRep:

```bash
m2c_hicrep="$m2c_core"
printf 'HICREP_ENV_PREFIX="%s"\nRUN_HICREP="true"\n' "$m2c_hicrep" >> "$m2c_setup_conf"
if bash scripts/bootstrap_environments.sh -c "$m2c_setup_conf" --check-only \
  2>&1 | tee "$m2c_work/audit/core-hicrep-check.log"; then
  printf '%s\n' 'Core also passed HiCRep; keep this prefix.'
else
  printf '%s\n' 'Review the log; create isolated HiCRep only if core passed and HiCRep failed.'
fi
```

Expected: both core and HiCRep probes succeed. The HiCRep probe checks its small self-concordance contract. If HiCRep alone is missing/incompatible, create a fresh owned HiCRep prefix:

```bash
m2c_hicrep="$m2c_base/envs/hicrep-0.1.0-shared-20261007"
printf 'HICREP_ENV_PREFIX="%s"\n' "$m2c_hicrep" >> "$m2c_setup_conf"
(umask 022; bash scripts/bootstrap_environments.sh -c "$m2c_setup_conf" --create-missing) \
  2>&1 | tee "$m2c_work/audit/hicrep-create.log"
bash scripts/bootstrap_environments.sh -c "$m2c_setup_conf" --check-only \
  2>&1 | tee "$m2c_work/audit/core-hicrep-recheck.log"
```

The separate environment uses `envs/hicrep.yml`; its pip package is installed only as part of that explicitly requested new environment creation. Do not add HiCRep to the upstream environment or the accepted core clone. Expected: both probes now succeed; existing environments remain unchanged.

## 9. Save the selected software inventory; defer differential R

```bash
printf 'role\tprefix\ncore\t%s\nhicrep\t%s\nbootstrap\t%s\nconda\t%s\n' \
  "$m2c_core" "$m2c_hicrep" "$m2c_bootstrap" "$m2c_conda" \
  > "$m2c_work/audit/chosen-prefixes.tsv"
"$m2c_conda" list -p "$m2c_core" --json > "$m2c_work/audit/core-packages.json"
"$m2c_conda" list -p "$m2c_hicrep" --json > "$m2c_work/audit/hicrep-packages.json"
cat "$m2c_work/audit/chosen-prefixes.tsv"
```

Confirm upstream Conda metadata and the selected home settings remain unchanged:

```bash
sha256sum --check "$m2c_work/audit/upstream-conda-meta-before.txt"
for m2c_home_file in "$HOME/.bashrc" "$HOME/.bash_profile" "$HOME/.profile" \
  "$HOME/.condarc" "$HOME/.conda/environments.txt"; do
  if [[ -f "$m2c_home_file" ]]; then
    sha256sum "$m2c_home_file"
  else
    printf 'ABSENT  %s\n' "$m2c_home_file"
  fi
done > "$m2c_work/audit/home-settings-after.txt"
cmp "$m2c_work/audit/home-settings-before.txt" "$m2c_work/audit/home-settings-after.txt"
```

If the source could not be inventoried in step 6, record that fact and omit its checksum command; do not fabricate a before-record. Differences require explanation before publication; do not overwrite another workflow's settings to make this check pass. Package metadata hashes detect recorded source changes but are not a complete byte-for-byte audit of every source environment file. Recheck home settings after the fixtures too.

Do not create the differential R environment for the initial descriptive installation. When inferential analysis is prepared, an existing ATAC/RNA **Conda** prefix may be reused only if its R/edgeR/limma/statmod/jsonlite checks, numerical tests and matching calibration pass. Otherwise use the new owned `envs/diff-r.yml` environment through explicit bootstrap. The current bootstrap does not manage an arbitrary system-R installation. Merely finding DESeq2 installed is insufficient for these contact/insulation models.

## 10. Run the small Linux smoke test, then the full descriptive fixture

Start with the six-sample QC/HiCRep/report fixture in a fresh directory:

```bash
m2c_smoke="$m2c_work/validation/qc-smoke-20261007"
"$m2c_conda" run --no-capture-output -p "$m2c_core" \
  python scripts/generate_demo.py --out "$m2c_smoke" --core-prefix "$m2c_core" --raw-only
printf 'CONDA_EXE="%s"\nPYTHON_BOOTSTRAP="%s"\nHICREP_ENV_PREFIX="%s"\n' \
  "$m2c_conda" "$m2c_bootstrap" "$m2c_hicrep" >> "$m2c_smoke/project.conf"
sha256sum "$m2c_smoke"/*.mcool > "$m2c_work/audit/smoke-input-sha256.txt"
bash scripts/preflight_check.sh -c "$m2c_smoke/project.conf" \
  2>&1 | tee "$m2c_work/audit/smoke-preflight.log"
bash scripts/microc2compare.sh -c "$m2c_smoke/project.conf" --dry-run
bash scripts/microc2compare.sh -c "$m2c_smoke/project.conf" \
  2>&1 | tee "$m2c_work/audit/smoke-run.log"
```

The generator's `--raw-only` disables cooltools features but retains **weighted QC**, HiCRep and reporting. It is not a no-weight mode. The demo preset reserves one worker/8,000 MB per task inside a 4-CPU/16,000-MB budget. If step 5 showed a smaller allocation, override global **and** task budgets in this generated configuration before preflight; see [resource guidance](07-server-resources.md).

Expected: six contact summaries, raw/balanced P(s) curves (contact frequency versus genomic distance), fifteen biological-sample HiCRep similarity comparisons and an HTML report:

```bash
m2c_smoke_run="$m2c_smoke/results/demo/synthetic_demo/run001"
ls "$m2c_smoke_run/08_reports/index.html"
bash scripts/microc2compare.sh -c "$m2c_smoke/project.conf" --resume
"$m2c_bootstrap" scripts/validate_outputs.py --config-json "$m2c_smoke_run/00_metadata/config.json"
sha256sum --check "$m2c_work/audit/smoke-input-sha256.txt"
```

Expected: resume succeeds, output validation succeeds, and every original map hash is `OK`. Review the report through your normal server-file browser or copy the complete `08_reports/` folder to your workstation; copying only the HTML loses relative figures/downloads.

Next verify the **full descriptive feature path**, using another fresh destination and omitting `--raw-only`:

```bash
m2c_full="$m2c_work/validation/full-descriptive-20261007"
"$m2c_conda" run --no-capture-output -p "$m2c_core" \
  python scripts/generate_demo.py --out "$m2c_full" --core-prefix "$m2c_core"
printf 'CONDA_EXE="%s"\nPYTHON_BOOTSTRAP="%s"\nHICREP_ENV_PREFIX="%s"\n' \
  "$m2c_conda" "$m2c_bootstrap" "$m2c_hicrep" >> "$m2c_full/project.conf"
sha256sum "$m2c_full"/*.mcool > "$m2c_work/audit/full-input-sha256.txt"
bash scripts/bootstrap_environments.sh -c "$m2c_full/project.conf" --check-only
bash scripts/preflight_check.sh -c "$m2c_full/project.conf"
bash scripts/microc2compare.sh -c "$m2c_full/project.conf" --dry-run
bash scripts/microc2compare.sh -c "$m2c_full/project.conf" \
  2>&1 | tee "$m2c_work/audit/full-descriptive-run.log"
bash scripts/microc2compare.sh -c "$m2c_full/project.conf" --resume
m2c_full_run="$m2c_full/results/demo/synthetic_demo/run001"
"$m2c_bootstrap" scripts/validate_outputs.py --config-json "$m2c_full_run/00_metadata/config.json"
sha256sum --check "$m2c_work/audit/full-input-sha256.txt"
```

Expected: enabled compartments, insulation/boundaries, loops/consensus, quantification, aggregates, descriptive comparisons and reports validate. Valid empty loop calls are possible; a failed caller is not a successful empty result. The fixture needs no R environment because its comparisons are descriptive. Fix any Linux/API/feature failure before a real-data run; a smoke-test pass alone does not validate the omitted features.

If you need smaller allocations, apply the resource overrides to **both** generated project files before their runs. Preserve logs and benchmark files. These tests establish installation behavior; production biological validity still requires the separate release gates.

## 11. Publish shared access, then inventory real microc2tracks inputs

Complete [shared publication and a non-maintainer smoke/full-fixture test](16-shared-server-installation.md) now. A maintainer-only pass is insufficient. That guide publishes a versioned software profile/command and explains ordinary-user setup without cloning code or creating duplicate environments. The remaining commands use this account's private workspace. Upstream data remain readable only by accounts authorized for those datasets; all-user software access does not require making private research data public.

Use the actual sheet and `OUTDIR` from a completed upstream run:

```bash
m2c_upstream_sheet="/ACTUAL/MICROC2TRACKS/config/samplesheet.csv"
m2c_upstream_results="/ACTUAL/MICROC2TRACKS/results"
"$m2c_conda" run --no-capture-output -p "$m2c_core" \
  python scripts/prepare_microc2compare_inputs.py \
  --samplesheet "$m2c_upstream_sheet" --upstream-results "$m2c_upstream_results" --preview
```

Expected: consolidated/single-source matrix suggestions, requested 10/100-kb groups, stored weights, related pairs/QC/browser/reference records and explicit issues. `needs_review` is not scientific acceptance. Missing resolution/weight/deduplication evidence requires upstream records or remediation, not downstream rebuilding.

Export one reference into a **new**, not-yet-created configuration directory:

```bash
m2c_project="$m2c_work/projects/first-real-pilot"
"$m2c_conda" run --no-capture-output -p "$m2c_core" \
  python scripts/prepare_microc2compare_inputs.py \
  --samplesheet "$m2c_upstream_sheet" --upstream-results "$m2c_upstream_results" \
  --project first_pilot --reference-id mm39 --write-drafts "$m2c_project/config"
```

Use your actual reference, e.g. `hg38` instead of `mm39`. Read `READ_BEFORE_RUNNING.md`. Complete donor/culture identity, technical-library relationships, original-source checksums, filtering/count units, balancing recipes and reference registration. Reuse existing assembly/chromosome/blacklist/phasing files only when their identities, coordinates and grids match; they are registered by path/checksum, not recopied by default. A phasing track is needed for compartments; do not invent one or silently orient compartments without it.

Follow [the handoff guide](13-upstream-handoff.md) to complete the drafts and configure a **descriptive** first pilot. Use the pilot preset only when the allocation permits it. Start with QC, inspect benchmarks, then resume appropriate descriptive modules. Pools, downsampling and new-map balancing remain off.

## 12. Keep the accepted shared deployment simple

The deployment has one central manager, one versioned source release, one accepted core clone/fallback, one HiCRep prefix only if needed, and one software profile/command. Every user supplies a private project configuration and uses these same checked paths. Keep accepted environments unchanged; provision/test a new version for an update. Use the scheduler or coordinate heavy jobs so this installation does not compete unexpectedly with your existing work.

| Item | Reuse when | Create separately when |
|---|---|---|
| Conda executable | Existing central installation outside homes works for ordinary accounts | Administrator provides a central manager only if none is suitable |
| Bootstrap Python | Existing Python is 3.11+ | No suitable existing interpreter exists |
| Core scientific software | Clone the existing upstream core once; accept only after checks/full Linux fixture | Use the fresh pinned core when cloning or compatibility fails |
| HiCRep prefix | Checked candidate passes HiCRep and fixture checks | HiCRep is absent/incompatible; do not add it to a shared prefix |
| Differential R prefix | Correct packages/numerics and matching calibration pass | Inferential analysis is needed and no suitable immutable prefix exists |
| Raw counts and weights | Existing biological-unit `.norm.mcool` is compatible | Correct missing/bad sample outputs upstream; default downstream never repairs them |
| References and tracks | Exact reference/axis/grid identities match | Supply the correct upstream/registered asset; only small analysis sidecars are newly written |
| CPU/RAM/storage | Your verified allocation and pilot measurements permit use | Adjust stage budgets/concurrency rather than reserving the entire machine |

Installation is ready for exploratory piloting when the selected prefixes, full Linux fixture, resume/output checks and unchanged map hashes pass and actual upstream registration is reviewed. Formal testing still requires its R/numerical/calibration/design/feature gates. Keep the [verification record](../validation/local-verification.md) distinction between already executed local checks and checks you still need to run on the server.

## If installation stops

| Failure | Next action |
|---|---|
| Conda/Python path missing | Confirm the site's installed executable/prefix; run the explicit version/path checks |
| Clone lacks a requested API/package | Correct only the documented missing-jsonschema case in the clone; otherwise use step 7, preserving the failed log |
| HiCRep missing | Create the isolated HiCRep prefix in step 8 |
| Solve/download fails | Preserve the log, check network/channel availability and solver information; do not alter shared packages or silently change pins |
| Permission error in a shared prefix/cache | Use an owned prefix/cache according to site policy; do not recursively change shared permissions |
| Disk/quota/allocation failure | Choose a suitable owned filesystem or permitted compute allocation; lower both global and task budgets as appropriate |
| Full descriptive fixture fails | Keep its logs/config/inventory and diagnose the specific adapter/runtime issue before real data |
| Old fixture/export directory exists | Choose a fresh destination; generators/export refuse overwriting |
| Two conditions lack independent units or have confounded metadata | Installation does not resolve experimental-design limitations; use an explicit descriptive path or correct the design |

The useful records to retain/share are `audit/server-inventory.txt`, `audit/chosen-prefixes.tsv`, package inventories, failed bootstrap logs and the relevant fixture task log/config. They describe the installed paths and failure without transferring large original maps. No server commands in this guide were executed during its preparation.
