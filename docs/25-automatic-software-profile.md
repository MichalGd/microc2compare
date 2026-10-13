# Automatic software settings for shared users

For a new installation from repository source, use [the GitHub installation guide](26-install-from-github.md). ZIP transfer sections here document earlier private deliveries; those artifacts are not included in Git.

Version 0.2.5 lets the maintainer bind a checked software profile to a versioned launcher. Once that launcher is activated, users run their own project directly:

```bash
microc2compare -c /home/USERNAME/my_project/project.conf --dry-run
microc2compare -c /home/USERNAME/my_project/project.conf
microc2compare -c /home/USERNAME/my_project/project.conf --analysis tads
```

No manual `source` line or Conda activation is required for these commands. Each project still chooses its input tables, `OUTDIR`, `SCRATCH_ROOT`, analysis parameters and resource budgets. Outputs are written under `OUTDIR/PROJECT/ANALYSIS_SET_ID/RUN_ID`; they do not go into the shared software directory. Users may select any directory they own, including their home. Since 0.3.1 runtime cache/temp paths are derived internally from `SCRATCH_ROOT/.runtime`; session exports are unnecessary for ordinary workflow commands. Since 0.3.2 [preparation and checking](32-simple-project-commands.md) also use the shared command.

## Which settings win?

The runner reads template defaults, then the trusted project configuration, then the selected software profile. Only these profile settings are imported:

```text
ENVIRONMENT_MODE, CONDA_EXE, PYTHON_BOOTSTRAP,
CORE_ENV_PREFIX, HICREP_ENV_PREFIX, DIFF_R_ENV_PREFIX, LOOP_ENV_PREFIX,
CORE_LOCK, HICREP_LOCK, DIFF_R_LOCK, LOOP_LOCK, R_SCRIPT
```

Their profile values override project software settings, including placeholders in a copied template. Profile assignments to input/output paths, resource budgets or analysis switches are not imported. Profiles are trusted Bash; this filtering is not a security sandbox.

A published launcher records the default profile's absolute path and SHA256. A missing or changed default profile stops the run before scientific analysis. Freeze the profile and its environment prefixes; publish a new pair when dependencies change. The normalized `00_metadata/config.json` records the effective software settings and `_software_profile` with its path, hash and selection (`launcher_default` or `explicit_override`).

The prepared server upgrade reuses `/opt/conda_envs/microc2compare/site/shared-r-ql-fix-20261008.conf`, the checked core/HiCRep prefixes and the existing `/opt/conda_envs/rna_ends2tracks-0.1.0a12.post1` R prefix. It installs no packages and changes no existing profile or environment.

## Explicit overrides and older runs

For an intentional software change, select a complete trusted profile explicitly:

```bash
microc2compare -c project.conf --software-profile /path/to/reviewed-software.conf --dry-run
```

Relative override paths resolve against the invocation directory. The equivalent environment variable is `M2C_SOFTWARE_PROFILE`; the command-line option takes precedence over it. An explicit override is hashed and recorded but does not need to match the launcher's default hash. Include an absolute `PYTHON_BOOTSTRAP` and all relevant software values in an override; omitted settings retain project/default values. A failed or unreadable override stops the command.

Existing project `source` lines continue to parse, but the selected launcher's software settings take precedence. Use the original versioned command and configuration for older runs. Direct `bash scripts/microc2compare.sh` invocation and unbound older launchers retain their existing explicit software configuration behavior; they can opt into the new profile option. Separate maintainer helpers still need their documented executable/profile settings. Automatic loading applies to the main launcher, its selected analyses and its internal stage helpers.

The shared configuration helper changed in this release, so conservative source fingerprints can invalidate previously completed stages. Resume may recompute them. Do not edit completion records to force reuse. Scientific algorithms and R adapters are unchanged from 0.2.4.

## Maintainer publication

Publish the new source/profile pair before activation:

```bash
m2c_base="/opt/conda_envs/microc2compare"
m2c_src="$m2c_base/releases/auto-profile-20261008/microc2compare-0.2.5"
m2c_site="$m2c_base/site/shared-r-ql-fix-20261008.conf"
/opt/miniconda/bin/python "$m2c_src/scripts/publish_shared_command.py" \
  --source "$m2c_src" --base "$m2c_base" --site-profile "$m2c_site"
```

The publisher verifies the release manifest, requires the profile below the shared `site` directory, rejects Unix group/world writable profiles or their immediate directory, and refuses to rebind an existing versioned launcher. Test `--list-analyses` under a non-owner account, then run a project without manual software lines before using the same publication command with `--activate`. `--help` alone does not load or validate the profile. Windows-local tests do not establish Linux permission/ACL acceptance.

## Upgrade this server

**On Windows PowerShell**, transfer the prepared delivery:

```powershell
$m2c_delivery = "D:\path\to\microc2compare\.local"
scp -P 22 "$m2c_delivery\microc2compare-0.2.5-auto-profile-8oct26.zip" "$m2c_delivery\microc2compare-0.2.5-auto-profile-8oct26.zip.sha256" "$m2c_delivery\install-microc2compare-0.2.5.py" "$m2c_delivery\install-auto-profile-8oct26.sh" "USERNAME@SERVER:/opt/conda_envs/microc2compare/incoming/"
```

**On server Bash**, verify the helper's SHA256 supplied with the delivery, then run:

```bash
bash /opt/conda_envs/microc2compare/incoming/install-auto-profile-8oct26.sh
```

The helper verifies and installs the new archive, checks that its R adapters match the already passed R audit, publishes the profile-bound versioned command, and tests profile loading under `nobody`. It creates a new private configuration beside the existing synthetic fixture with the manual software lines removed, runs QC through the new launcher, validates outputs/profile provenance and checks the original matrix hashes. It activates the convenience command only after those checks succeed. Existing projects, releases, profiles and environments stay intact; retries accept only matching existing delivery/configuration content.

Expected final message:

```text
Automatic software profile accepted and activated.
Users can run: microc2compare -c project.conf
```

This interface check uses six synthetic `chrDemo` maps, not real mouse/human data or statistical calibration. The maintainer subsequently reported valid QC/report outputs, automatic-profile provenance verification, unchanged matrices/original configuration and shared 0.2.5 activation. Ordinary-researcher full execution, intended-design calibration and a real-data pilot remain separate validation steps. See [the verification record](../validation/local-verification.md).
