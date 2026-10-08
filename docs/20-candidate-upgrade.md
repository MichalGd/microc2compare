# Install and accept the modules/TAD candidate

For a new installation from repository source, use [the GitHub installation guide](26-install-from-github.md). ZIP transfer sections here document earlier private deliveries; those artifacts are not included in Git.

Version 0.2.0 has a bootstrap-state import error discovered during server acceptance. Use [the 0.2.1 correction](21-bootstrap-state-correction.md) before accepting/activating this candidate. These original 0.2.0 transfer/install steps remain the record for the installed candidate; existing demo inputs can be resumed without regeneration.

This installs version 0.2.0 beside the accepted 0.1.0 release. It reuses the checked shared core clone and HiCRep prefix. No additional package is required for the insulation-based TAD method. GitHub creation is deferred. The R environment and inferential calibration are separate optional work described in [the R guide](19-r-environment-and-inference.md).

The new source must pass Linux scientific/module acceptance before activating the shared alias. Previous owner-account success on 0.1.0 is not proof that this candidate passed. Results, audit files, scratch and caches belong to each user's chosen directory; there is no required shared analysis/results folder.

## 1. Transfer from Windows

Run these in **Windows PowerShell**, where the prompt begins `PS`, not `USERNAME@SERVER`. The prepared local delivery includes the archive, checksum sidecar and standalone installer. Assign the local directory first:

```powershell
$m2c_delivery = "D:\path\to\microc2compare\.local"
scp -P 22 "$m2c_delivery\microc2compare-0.2.0-tads-modules-8oct26.zip" "$m2c_delivery\microc2compare-0.2.0-tads-modules-8oct26.zip.sha256" "$m2c_delivery\install-microc2compare-0.2.0.py" "USERNAME@SERVER:/opt/conda_envs/microc2compare/incoming/"
```

Use the same authentication as the successful earlier transfers. Lowercase `-p` is not the SSH port flag. Wait for a successful transfer before the server commands.

## 2. Start a dedicated Linux session

Run **this command alone in the server terminal**:

```bash
bash --noprofile --norc
```

Then restore the dedicated-session variables. This affects only this shell; no startup file is edited. A personal project directory other than the suggested HOME workspace is allowed.

```bash
unset HISTFILE BASH_ENV ENV PYTHONHOME PYTHONPATH R_HOME R_LIBS R_LIBS_USER
umask 077
m2c_base="/opt/conda_envs/microc2compare"
m2c_work="$HOME/microc2compare-work"
m2c_conda="/opt/miniconda/bin/conda"
m2c_bootstrap="/opt/miniconda/bin/python"
mkdir -p "$m2c_work"/{audit,validation,tmp,runtime-config,runtime-data}
mkdir -p "$m2c_work/cache"/{conda-pkgs,matplotlib,pip}
export TMPDIR="$m2c_work/tmp"
export CONDA_PKGS_DIRS="$m2c_work/cache/conda-pkgs"
export CONDA_REGISTER_ENVS=false CONDA_ALWAYS_COPY=true
export XDG_CACHE_HOME="$m2c_work/cache"
export XDG_CONFIG_HOME="$m2c_work/runtime-config"
export XDG_DATA_HOME="$m2c_work/runtime-data"
export MPLCONFIGDIR="$m2c_work/cache/matplotlib"
export PIP_CACHE_DIR="$m2c_work/cache/pip"
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
```

## 3. Verify and extract a new source release

These are **Linux Bash commands**, one step at a time. The installer verifies the whole archive and every manifest-listed file before writing the new release. It refuses an existing release ID and does not touch an environment, project or command link.

```bash
cd "$m2c_base/incoming"
sha256sum --check microc2compare-0.2.0-tads-modules-8oct26.zip.sha256
```

Continue only if the checksum is `OK`:

```bash
read -r m2c_sha m2c_archive_name < microc2compare-0.2.0-tads-modules-8oct26.zip.sha256
"$m2c_bootstrap" "$m2c_base/incoming/install-microc2compare-0.2.0.py" \
  --archive "$m2c_base/incoming/microc2compare-0.2.0-tads-modules-8oct26.zip" \
  --sha256 "$m2c_sha" --base "$m2c_base" --release-id tads-modules-20261008
m2c_src="$m2c_base/releases/tads-modules-20261008/microc2compare-0.2.0"
bash "$m2c_src/scripts/microc2compare.sh" --help
PYTHON_BOOTSTRAP="$m2c_bootstrap" bash "$m2c_src/scripts/microc2compare.sh" --list-analyses
```

Expected: a verified new source path, help including `--analysis`, and the catalogue. If anything fails, stop at that step and inspect the error. Do not remove the accepted release to make room for a failed candidate.

## 4. Run acceptance in your own workspace

The existing software-only profile points to the already checked core/HiCRep environments and leaves R unconfigured. Reusing that profile does not force the old workflow source. It contains no dataset/output path.

```bash
m2c_site="$m2c_base/site/shared-design-json-fix-20261008.conf"
m2c_accept="$m2c_work/validation/modules-candidate-20261008"
set -o pipefail
bash "$m2c_src/scripts/validate_module_demo.sh" \
  --site "$m2c_site" --out "$m2c_accept" \
  2>&1 | tee "$m2c_work/audit/modules-candidate-20261008.log"
m2c_accept_status=${PIPESTATUS[0]}
printf 'Candidate acceptance exit code: %s\n' "$m2c_accept_status"
```

The helper refuses an existing destination, generates original artificial data, and executes the full descriptive pipeline, resume, eleven named selections and a boundary-threshold variant. It validates outputs and compares matching compartment/insulation/loop tables between full and selected executions. Finally it verifies original matrix hashes and that results contain no new `.cool`/`.mcool`. It uses the four-CPU/16,000-decimal-MB demo preset; confirm that allocation first. On a smaller allocation, supply `--resource-config /owned/lower-resources.conf` with lower global **and** stage reservations rather than editing accepted software. This is a trusted Bash resource preset and is applied before the software-only profile.

Expected final line: `Full/selected/variant acceptance passed`. An empty TAD/event family is a valid software result; known-event tests and biological inspection are separate. Inspect the full HTML plus selected TAD/variant reports, coverage flags and downloadable tables. Optional pool/depth creation is not part of this default acceptance, and formal R/calibration is not claimed.

## 5. Repeat as an ordinary researcher

The second account uses its own dedicated session and its own `m2c_work`. That account must be able to read the shared profile/source and execute the checked prefixes, while lacking write access to shared software. Run the same acceptance helper with a new destination under that account's own chosen workspace. Do not make USERNAME's private inputs or results public to perform this test. See [ordinary-account checks](16-shared-server-installation.md#4-accept-the-installation-from-a-second-account).

## 6. Publish the accepted command

Only after acceptance and report review, the maintainer publishes a versioned launcher. This first command does not change the shared alias:

```bash
"$m2c_bootstrap" "$m2c_src/scripts/publish_shared_command.py" --source "$m2c_src" --base "$m2c_base"
"$m2c_base/bin/microc2compare-0.2.0-tads-modules-20261008" --help
```

Activate the accepted candidate explicitly:

```bash
"$m2c_bootstrap" "$m2c_src/scripts/publish_shared_command.py" --source "$m2c_src" --base "$m2c_base" --activate
hash -r
microc2compare --help
```

The existing `/usr/local/bin/microc2compare` link continues to point to the shared alias, so no additional sudo operation is needed. Old versioned launchers/releases remain available for reproducibility. Use the old pinned command for old runs; use new run IDs for the new release. Source changes can legitimately invalidate old completion signatures.

## 7. Download and open the report on Windows

Run these in **Windows PowerShell**, one line at a time. The explicit paths avoid dependence on variables from an earlier window:

```powershell
New-Item -ItemType Directory -Force -Path "$env:USERPROFILE\Downloads\microc2compare-0.2.0-demo" | Out-Null
scp -P 22 -r "USERNAME@SERVER:/home/USERNAME/microc2compare-work/validation/modules-candidate-20261008/results/demo/synthetic_demo/run001" "$env:USERPROFILE\Downloads\microc2compare-0.2.0-demo"
Start-Process -FilePath "msedge.exe" -ArgumentList ('"' + $env:USERPROFILE + '\Downloads\microc2compare-0.2.0-demo\run001\08_reports\index.html"')
```

Opening a browser normally prints no terminal output. You can also press Ctrl+O in the browser and select the downloaded `index.html`. Download the entire run folder so figures and result links remain usable.

## Remaining scientific work

Pilot actual assembly-matched mouse/human data and inspect reference assets, replicate identity, masking, domain sensitivity, comparisons and resource measurements. Finish R numerical checks/calibration only when inference is needed. Candidate installation, scientific validation and GitHub publication are separate milestones. Record actual outcomes in your private audit, including failed checks; neither a local archive nor a successful synthetic run establishes all biological claims.
