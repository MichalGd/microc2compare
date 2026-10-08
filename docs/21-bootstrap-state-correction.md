# Correct bootstrap state checks and resume acceptance

For a new installation from repository source, use [the GitHub installation guide](26-install-from-github.md). ZIP transfer sections here document earlier private deliveries; those artifacts are not included in Git.

**Latest combined candidate:** 0.2.2 includes this correction plus differential loop volcano plots. Use [the amended installation guide](22-loop-volcano-plots.md) directly; the 0.2.1 instructions below retain the correction history.

The 0.2.0 server acceptance stopped because registration-signature calculation imported `matrix.py` under the bootstrap Python, which intentionally does not need scientific packages. Version 0.2.1 moves the unchanged grid-selection policy into the dependency-free analysis planner. The matrix adapters and state checks share that policy. Existing core/HiCRep environments and stored upstream weights need no changes.

The new source is installed beside both previous releases. Do not overwrite their files/manifests. The global command remains on the previously accepted release until the new acceptance checks pass. GitHub and formal R/calibration work remain separate.

## 1. Transfer from Windows PowerShell

Run where the prompt contains `PS`, not in the server Bash terminal:

```powershell
$m2c_delivery = "D:\path\to\microc2compare\.local"
scp -P 22 "$m2c_delivery\microc2compare-0.2.1-bootstrap-fix-8oct26.zip" "$m2c_delivery\microc2compare-0.2.1-bootstrap-fix-8oct26.zip.sha256" "$m2c_delivery\install-microc2compare-0.2.1.py" "USERNAME@SERVER:/opt/conda_envs/microc2compare/incoming/"
```

## 2. Verify and install in server Bash

Use the dedicated session and private cache/temp settings from [candidate installation](20-candidate-upgrade.md#2-start-a-dedicated-linux-session). Restore the following variables in that session:

```bash
m2c_base="/opt/conda_envs/microc2compare"
m2c_work="$HOME/microc2compare-work"
m2c_bootstrap="/opt/miniconda/bin/python"
cd "$m2c_base/incoming"
sha256sum --check microc2compare-0.2.1-bootstrap-fix-8oct26.zip.sha256
```

Continue only if the checksum reports `OK`:

```bash
read -r m2c_sha m2c_archive_name < "$m2c_base/incoming/microc2compare-0.2.1-bootstrap-fix-8oct26.zip.sha256"
"$m2c_bootstrap" "$m2c_base/incoming/install-microc2compare-0.2.1.py" \
  --archive "$m2c_base/incoming/microc2compare-0.2.1-bootstrap-fix-8oct26.zip" \
  --sha256 "$m2c_sha" --base "$m2c_base" --release-id bootstrap-fix-20261008
```

The installer verifies every manifest entry before writing, refuses an existing release ID and changes neither packages nor the shared command. After successful installation:

```bash
m2c_src="$m2c_base/releases/bootstrap-fix-20261008/microc2compare-0.2.1"
bash "$m2c_src/scripts/microc2compare.sh" --help
bash "$m2c_src/scripts/validate_module_demo.sh" --help
```

## 3. Resume the existing acceptance fixture

The helper's new `--resume` verifies the original checksum list before running anything scientific. It requires the original `project.conf`, checksum file and the same software-profile path. It preserves existing configuration files and creates missing standalone configurations only. It does not invoke the demo generator or reset input hashes. Omit a new `--resource-config` when resuming; the existing project retains its resource settings.

```bash
m2c_site="$m2c_base/site/shared-design-json-fix-20261008.conf"
m2c_accept="$m2c_work/validation/modules-candidate-20261008"
set -o pipefail
bash "$m2c_src/scripts/validate_module_demo.sh" \
  --site "$m2c_site" --out "$m2c_accept" --resume \
  2>&1 | tee "$m2c_work/audit/modules-candidate-bootstrap-fix-20261008.log"
m2c_accept_status=${PIPESTATUS[0]}
printf 'Candidate acceptance exit code: %s\n' "$m2c_accept_status"
```

Changed source hashes can legitimately recompute earlier tasks. Success requires `Full/selected/variant acceptance passed` and exit code zero, valid output records, agreement checks and unchanged original matrices. This is actual server acceptance, not something established by packaging or local tests. If a task fails, retain its log and the existing fixture; no hand editing of completion signatures or blanket lock removal is appropriate.

For a researcher who has no existing fixture, omit `--resume` and supply a new private destination. Repeat from the ordinary account's own workspace, review reports and follow [publication after acceptance](20-candidate-upgrade.md#6-publish-the-accepted-command), using the 0.2.1 `m2c_src`. Its versioned launcher is `microc2compare-0.2.1-bootstrap-fix-20261008`; retain the previous launcher/releases.
