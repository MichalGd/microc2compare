# Install from GitHub on a shared Linux server

This guide starts from repository source. Historical upgrade guides describe previously prepared local ZIP deliveries; those private artifacts are not included in Git. If the shared command is already installed, ordinary users can skip installation and go to step 5.

One maintainer manages software outside personal homes. Each researcher controls their own inputs, configuration, results and scratch. Analyses never install or update packages.

## 1. Confirm existing shared software

Ask the administrator for a writable shared application directory and a central Conda manager. Examples use `/opt/conda_envs/microc2compare` and `/opt/miniconda`. Do not change permissions on upstream installations or install shared software in a private home.

Use Linux Bash 4.4+, bootstrap Python 3.11+, Git and `sha256sum`. Inspect and test existing core/HiCRep environments before provisioning anything. A checked independent clone of microc2tracks core is the recommended reuse path. A checked shared R/edgeR/limma prefix can be reused unchanged; R is optional for descriptive analysis. See [detailed environment preparation](15-server-installation-start.md), [environment policy](01-installation.md) and [R numerical checks](19-r-environment-and-inference.md).

## 2. Clone and install an immutable source release

**Maintainer, Linux Bash:** adjust base/bootstrap paths to your site. For a private repository, authenticate through normal GitHub/Git credentials first. Select a reviewed commit/tag when available; `main` can change. Choose another new checkout destination if the example already exists.

```bash
set -Eeuo pipefail
umask 022
m2c_base="/opt/conda_envs/microc2compare"
m2c_bootstrap="/opt/miniconda/bin/python"
mkdir -p "$m2c_base"/{incoming,releases,envs,site,bin}
m2c_checkout="$m2c_base/incoming/github-source-20261009"
test ! -e "$m2c_checkout"
git clone https://github.com/MichalGd/microc2compare.git "$m2c_checkout"
git -C "$m2c_checkout" log -1 --oneline
```

Review the commit and validation status before continuing. The archive captures this checkout; never `git pull` inside a deployed release.

```bash
m2c_revision=$(git -C "$m2c_checkout" rev-parse HEAD)
m2c_release="github-${m2c_revision:0:12}"
m2c_version=$("$m2c_bootstrap" -c 'import sys,tomllib; print(tomllib.load(open(sys.argv[1],"rb"))["project"]["version"])' "$m2c_checkout/pyproject.toml")
m2c_archive="$m2c_base/incoming/microc2compare-$m2c_version-$m2c_release.zip"
"$m2c_bootstrap" "$m2c_checkout/scripts/package_release.py" --out "$m2c_archive"
sha256sum "$m2c_archive" > "$m2c_archive.sha256"
read -r m2c_sha m2c_archive_name < "$m2c_archive.sha256"
"$m2c_bootstrap" "$m2c_checkout/scripts/install_release.py" \
  --archive "$m2c_archive" --sha256 "$m2c_sha" \
  --base "$m2c_base" --release-id "$m2c_release"
m2c_src="$m2c_base/releases/$m2c_release/microc2compare-$m2c_version"
```

Expected: all allowlisted file hashes verify and the installer prints the source path. Existing archive/release destinations are refused. Preserve older versions for active runs.

## 3. Select and check environments once

Create a versioned software profile from the template. Replace **every** `/ACTUAL/...` value with a compatible existing executable/prefix. Set `DIFF_R_ENV_PREFIX=""` if installing descriptive analyses only.

```bash
m2c_site="$m2c_base/site/software-$m2c_release.conf"
test ! -e "$m2c_site"
cp "$m2c_src/config/site.example.conf" "$m2c_site"
chmod 0644 "$m2c_site"
```

Edit the profile before continuing. Keep only software settings. The profile, source and needed environment files must be readable/executable by intended accounts and not writable by ordinary users. Preserve Conda binary permissions; do not recursively apply text-file modes to environments. Freeze the profile after publication.

**Maintainer checks use a private workspace:**

```bash
umask 077
m2c_check="$HOME/microc2compare-work/validation/$m2c_release"
mkdir -p "$m2c_check"
test ! -e "$m2c_check/setup.conf"
printf 'source "%s"\nRUN_DIFFERENTIAL="false"\n' "$m2c_site" > "$m2c_check/setup.conf"
bash "$m2c_src/scripts/bootstrap_environments.sh" -c "$m2c_check/setup.conf" --check-only
bash "$m2c_src/scripts/validate_module_demo.sh" \
  --site "$m2c_site" --out "$m2c_check/demo"
```

Expected: package probes and full/resume/selected/TAD-variant descriptive acceptance pass, with unchanged matrix hashes. The fixture defaults to four CPUs/16,000 decimal MB; use `--resource-config` with smaller global **and task** budgets when needed. `--check-only` never installs packages. Failed prefixes need separate reviewed replacements. Optional R needs [actual adapter checks](19-r-environment-and-inference.md), not just package imports.

## 4. Publish, check access, then activate

```bash
"$m2c_bootstrap" "$m2c_src/scripts/publish_shared_command.py" \
  --source "$m2c_src" --base "$m2c_base" --site-profile "$m2c_site"
m2c_command="$m2c_base/bin/microc2compare-$m2c_version-$m2c_release"
"$m2c_command" --list-analyses
"$m2c_command" -c "$m2c_check/demo/project.conf" --dry-run
```

The launcher pins the profile's path/SHA256. Test `--list-analyses` under another account and run a private project through it without manual software lines before activation. `--help` alone does not test profile loading. Follow [ordinary-account permission/full-run checks](16-shared-server-installation.md#4-accept-the-installation-from-a-second-account); root/owner execution is not ordinary-account acceptance.

```bash
"$m2c_bootstrap" "$m2c_src/scripts/publish_shared_command.py" \
  --source "$m2c_src" --base "$m2c_base" --site-profile "$m2c_site" --activate
```

An existing `/usr/local/bin/microc2compare` link to `$m2c_base/bin/microc2compare` needs no replacement. Otherwise ask the administrator to create that link once after confirming the destination is unused, or invoke the absolute versioned command. Preserve old releases/profiles/prefixes.

## 5. Start a private real-data project

Use the [input helper](13-upstream-handoff.md) to export and review drafts. Confirm biological/technical identity, reference assembly and processing evidence. Start with explicit descriptive comparisons. Choose `OUTDIR` and `SCRATCH_ROOT` in your own writable directories.

```bash
microc2compare -c /path/to/your/project.conf --dry-run
microc2compare -c /path/to/your/project.conf --until-stage qc
# After inspecting QC and resource use:
microc2compare -c /path/to/your/project.conf --resume
# Or select an analysis in an appropriately configured run:
microc2compare -c /path/to/your/tad-project.conf --analysis tads
```

No software-profile line or Conda activation is needed for the bound command. Private cache/temp exports are optional [session setup](16-shared-server-installation.md#3-start-as-an-ordinary-server-user). The [module catalogue](17-analysis-modules.md) lists inputs, parameters and outputs; open `08_reports/index.html`.

Statistical calibration is separate from installation and matrix balancing: it tests false discoveries on simulated datasets with known truth. Descriptive analyses do not need it. Formal/preliminary differential loops/insulation need matching reviewed calibration, eligible biological designs and the required feature universe. Read [statistics](05-statistics.md) and [the protocol](../validation/calibration-protocol.md) before enabling inference.
