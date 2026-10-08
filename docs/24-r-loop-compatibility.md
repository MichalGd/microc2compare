# R loop adapter correction and server acceptance

For a new installation from repository source, use [the GitHub installation guide](26-install-from-github.md). ZIP transfer sections here document earlier private deliveries; those artifacts are not included in Git.

Version **0.2.4** fixes the loop adapter failure observed when testing the existing shared RNA environment. The active descriptive 0.2.3 installation and its accepted demonstration remain available. Install this candidate beside it; test R before activating a new shared command.

## Cause and scope

The user reported successful loading of edgeR 4.4.0, limma 3.62.1, statmod 1.5.2 and jsonlite 2.0.0 from `/opt/conda_envs/rna_ends2tracks-0.1.0a12.post1`. The loop adapter then failed with `vapply ... result is length 0` while computing approximate standard errors.

edgeR 4.4 stores moderated quasi-likelihood dispersions in `s2.post`; the adapter read the older `var.post` field, which was absent. The corrected helper uses `s2.post` when present and supports the older field name, while rejecting absent, nonfinite or malformed values. It also expands a common scalar negative-binomial dispersion across features and uses edgeR's actual working dispersion (`dispersion / average.ql.dispersion` for the adjusted-deviance fit) for information-based uncertainty and NB Pearson residuals. Unsupported dispersion matrices fail explicitly.

The edgeR fitting/testing calls, count filtering, exposure offsets, log-fold changes, P values and BH adjustment are unchanged. The approximate confidence intervals remain pointwise approximations requiring coverage validation. This correction does not enable inference or provide calibration evidence. See the [edgeR manual](https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf) and the upstream [QL implementation](https://raw.githubusercontent.com/bioc/edgeR/RELEASE_3_20/R/glmQLFTest.R) and [GLM implementation](https://raw.githubusercontent.com/bioc/edgeR/RELEASE_3_20/R/glmfit.R).

New loop-result columns record `working_NB_dispersion`, `QL_dispersion` and `QL_dispersion_field`; `dispersion` retains the dispersion stored by edgeR. The helper is included in the differential source fingerprint and R audit evidence. Earlier differential completion/calibration records are not compatible with this changed source. Other descriptive scientific stage source fingerprints are unchanged from 0.2.3; normal prerequisite checks still govern reuse.

Local numerical tests use an isolated project-owned Windows runtime with R 4.4.3, edgeR 4.4.2 and limma 3.62.2. They compare the adapter with direct edgeR results and an independent analytical two-group variance, check exposure scaling, and exercise the real insulation adapter. Base-R tests cover both field names, common/per-feature dispersion, scaling and invalid shapes/values. These local versions are distinct from the server candidate: Linux acceptance with edgeR 4.4.0 remains required.

## 1. Transfer from Windows PowerShell

Run at the Windows **PS** prompt:

```powershell
$m2c_delivery = "D:\path\to\microc2compare\.local"
scp -P 22 "$m2c_delivery\microc2compare-0.2.4-r-ql-fix-8oct26.zip" "$m2c_delivery\microc2compare-0.2.4-r-ql-fix-8oct26.zip.sha256" "$m2c_delivery\install-microc2compare-0.2.4.py" "USERNAME@SERVER:/opt/conda_envs/microc2compare/incoming/"
```

## 2. Verify and install in server Bash

Use the dedicated Bash session with private caches/temp from the [installation guide](20-candidate-upgrade.md#2-start-a-dedicated-linux-session). Run at the server **bash-5.2$** prompt:

```bash
m2c_base="/opt/conda_envs/microc2compare"
m2c_bootstrap="/opt/miniconda/bin/python"
cd "$m2c_base/incoming"
sha256sum --check microc2compare-0.2.4-r-ql-fix-8oct26.zip.sha256
```

Continue only after checksum `OK`:

```bash
read -r m2c_sha m2c_archive_name < "$m2c_base/incoming/microc2compare-0.2.4-r-ql-fix-8oct26.zip.sha256"
"$m2c_bootstrap" "$m2c_base/incoming/install-microc2compare-0.2.4.py" \
  --archive "$m2c_base/incoming/microc2compare-0.2.4-r-ql-fix-8oct26.zip" \
  --sha256 "$m2c_sha" --base "$m2c_base" --release-id r-ql-fix-20261008
```

The installer refuses an existing release ID and verifies all source hashes. It changes no old release, environment, project or command alias. After installation:

```bash
m2c_src="$m2c_base/releases/r-ql-fix-20261008/microc2compare-0.2.4"
bash "$m2c_src/scripts/microc2compare.sh" --help
```

## 3. Repeat the actual R numerical checks

Use a new private audit directory; retain the failed 0.2.3 log for comparison. No packages are added to the RNA environment.

```bash
umask 077
m2c_conda="/opt/miniconda/bin/conda"
m2c_core="/opt/conda_envs/microc2compare/envs/core-0.1.0-clone-20261007"
m2c_r_candidate="/opt/conda_envs/rna_ends2tracks-0.1.0a12.post1"
m2c_r_audit="$HOME/microc2compare-work/validation/r-024-check-$(date +%Y%m%d-%H%M%S)"
unset R_HOME R_LIBS R_LIBS_USER
export CONDA_REGISTER_ENVS=false

"$m2c_conda" run --no-capture-output -p "$m2c_core" \
  python "$m2c_src/scripts/validate_r_environment.py" \
  --conda "$m2c_conda" --prefix "$m2c_r_candidate" --out "$m2c_r_audit"
m2c_r_status=$?
printf '\nR validation exit code: %s\nAudit directory: %s\n' \
  "$m2c_r_status" "$m2c_r_audit"
```

Expected: exit zero, `R package and numerical checks passed; statistical calibration remains required`, and `acceptance.json` containing `passed_numerical_checks_only` / `not_run`. If it fails, inspect the named audit log. Do not modify the upstream environment, replace missing variances with one, skip uncertainty checks, or edit acceptance status.

## 4. Software profile, calibration and publication

Only after numerical acceptance, inspect ordinary-user access and the candidate prefix's immutability before choosing direct reuse or an isolated clone. Keep the existing descriptive profile; record the accepted R prefix in a new versioned software-only profile. Each researcher chooses their own project input/output/scratch directories. Follow [R environment and statistical acceptance](19-r-environment-and-inference.md) for the intended comparison design, calibration and biological pilot. Local tests establish neither calibration nor compatibility of arbitrary R package versions.

The maintainer can publish a versioned command after candidate acceptance, then explicitly activate it once the deployment checks are complete, using [shared command publication](23-compartment-agreement.md#4-publish-after-acceptance-and-ordinary-account-checks). Leave the working 0.2.3 alias active while testing this candidate. No GitHub creation or publication occurs during installation.
