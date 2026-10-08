# Shared installation for all server users

**Current shared-command usage:** a profile-bound 0.2.5 launcher loads software settings internally. Use [the automatic-profile upgrade guide](25-automatic-software-profile.md) for this server. The detailed manual profile/publication commands below describe direct helper invocation and older unbound deployments.

The supported server layout is **one shared software installation, separate user workspaces**. The runner already resolves its own source directory and writes run locks, logs, outputs and scratch outside it when configured that way. The recommended core is a one-time versioned clone of microc2tracks, with changes confined to that clone; an existing central Conda manager is reused. Your home installation and upstream workflows stay in place. Python bytecode writes are disabled by the Bash helpers so ordinary runs do not populate a shared release or prefix with caches.

All server accounts may read/run published code and checked environments. Only the designated installation maintainer may change them. Each account controls its own project configuration, manifests, results, temporary files and caches. Permission to use the workflow is separate from permission to read a particular research dataset.

Run the role-specific steps in order and stop at a failed command. Do not bypass a destination-exists guard: inspect an existing installation or select a new versioned destination instead. Commands below are instructions for Linux deployment; they have not been executed remotely.

## Layout and responsibilities

Example paths below assume `/opt/conda_envs/microc2compare`; substitute the approved shared application filesystem consistently. Do not choose a personal home for shared software.

```text
/opt/conda_envs/microc2compare/                         maintainer-owned; all users can traverse
  releases/shared-design-json-fix-20261008/microc2compare-0.1.0/
    scripts/, src/, config/, docs/, ...     published source; users read/run only
  envs/                                    one core clone; isolated HiCRep only if needed
  site/shared-design-json-fix-20261008.conf                 checked software paths; no project data
  bin/microc2compare-0.1.0-shared-design-json-fix-20261008    command pinned to this release
  bin/microc2compare                        convenience command for new runs

$HOME/microc2compare-work/                  private to each account
  projects/PROJECT/config/                  trusted .conf and biological tables
  projects/PROJECT/results/                 outputs, logs, reports and run locks
  projects/PROJECT/scratch/                 task-owned temporary scientific files
  validation/, audit/, tmp/, cache/         this user's checks and runtime files

/approved/upstream/results/                 existing maps; authorized read-only access
/approved/shared/references/                matching registered reference assets
```

Each user works in their own home directory, by default `$HOME/microc2compare-work`. For USERNAME this normally resolves to `/home/USERNAME/microc2compare-work`; other users get their own home paths. A different personal project directory is optional. Only workflow code, shared environments and site launchers/profiles belong outside home; `/dysk3` is not required. A collaborative project can use a project-specific group/ACL, but one account should operate each run. Never use a common writable `config/project.conf`, common demo destination or shared scratch directory for unrelated users. The existing run lock rejects simultaneous writers to the same run; it is not a server-wide job queue.

Conda documents centrally managed, read-only installations and user fallback locations in its [multi-user administration guide](https://docs.conda.io/projects/conda/en/stable/user-guide/configuration/admin-multi-user-install.html). Reuse an existing site installation where it passes the checks below; do not install a second Conda manager by default.

## 1. Provision the shared software root once

**Administrator, only if the root does not already exist:** replace `MAINTAINER` and `MAINTAINER_GROUP` with the actual installation account/group. The following creates the named directory; it does not change ownership of an existing upstream installation.

```bash
sudo install -d -o MAINTAINER -g MAINTAINER_GROUP -m 0755 /opt/conda_envs/microc2compare
```

Ensure every parent permits directory traversal by intended accounts. Confirm site ACLs and filesystem mounts permit reading files and running the selected Conda/Python binaries. Software may be readable by all accounts without making project data readable by them. If the software root already exists, review its owner/mode/ACL before changing anything; do not recursively change `/srv`, upstream environments, package caches or home directories.

**Maintainer:** follow [the numbered installation guide](15-server-installation-start.md), steps 1–10. It keeps shared software under `$m2c_base`, private audits/tests under `$m2c_work`, and uses `umask 022` only for creating shared software. Private work uses `umask 077`.

Use a central Conda executable and bootstrap Python outside personal homes, including their resolved targets. The upstream core is copied once with `conda create --clone --copy` into its final shared destination; ordinary filesystem copying/moving is not the method. Review copied activation hooks, editable packages and external symlinks for home dependencies. If cloning is unsuitable, use the clean fallback. Never move your existing home installation, alter its packages/startup settings or broaden its permissions. If no central manager is available, the administrator provides an independent one; the private installation remains untouched.

## 2. Publish checked software paths and a command

**Maintainer, after the installation checks:** keep the variables selected in the numbered guide. The following sanity checks catch paths missing from this account; second-account checks follow separately.

```bash
test -r "$m2c_src/scripts/microc2compare.sh"
test -x "$m2c_conda"
test -x "$m2c_bootstrap"
test -d "$m2c_core/conda-meta"
test -d "$m2c_hicrep/conda-meta"
```

For the newly extracted, owned source tree only, make directories traversable and source files readable by everyone, without granting group/other write access. This archive contains source text, not executable scientific binaries; the scripts are invoked with Bash/Python. Preserve executable permissions of Conda binaries rather than applying these source-only modes to environments.

```bash
test "$m2c_src" = "$m2c_base/releases/shared-design-json-fix-20261008/microc2compare-0.1.0" && \
  find "$m2c_src" -type d -exec chmod 0755 {} + && \
  find "$m2c_src" -type f -exec chmod 0644 {} +
```

Check the owner/group and access of each selected software prefix and its parent directories. Fresh prefixes were created under the shared-software umask; site ACLs can still affect their accessibility. Do not repair an inaccessible reused prefix by recursive chmod/chown. Select an approved shared prefix or provision a fresh one instead. Inspect activation hooks too: they must not require writes to shared software paths. The non-maintainer test is the acceptance check.

For the **new application-owned environments only**, establish all-user read/execute access and remove group/other write access. These prefixes were created with package copies. The guard resolves the target and confines changes to this application's environment directory; an existing external bootstrap/prefix retains its original permissions. Review ACLs and external link targets separately:

```bash
m2c_envs_root=$(readlink -f "$m2c_base/envs")
for m2c_prefix in "$m2c_core" "$m2c_hicrep" "$(dirname "$(dirname "$m2c_bootstrap")")"; do
  m2c_resolved_prefix=$(readlink -f "$m2c_prefix")
  case "$m2c_resolved_prefix" in
    "$m2c_envs_root"/*)
      test ! -L "$m2c_prefix"
      test -O "$m2c_resolved_prefix"
      chmod -R a+rX,go-w -- "$m2c_resolved_prefix"
      ;;
    *) printf 'Preserving existing external prefix permissions: %s\n' "$m2c_prefix" ;;
  esac
done
```

Do this during initial provisioning, before publication; never apply it to an environment already serving active runs. Do not extend the guard to `/opt/conda_envs`, upstream prefixes, Conda base or home directories. This establishes software permissions only and does not grant access to research inputs.

Create a versioned software profile outside the release. [The site template](../config/site.example.conf) documents these fields. Leave R unconfigured during descriptive setup; later inference uses a new reviewed profile with its checked R prefix and calibration evidence.

```bash
m2c_site="$m2c_base/site/shared-design-json-fix-20261008.conf"
test ! -e "$m2c_site"
test ! -L "$m2c_site"
cat > "$m2c_site" <<EOF
ENVIRONMENT_MODE="reuse"
CONDA_EXE="$m2c_conda"
PYTHON_BOOTSTRAP="$m2c_bootstrap"
CORE_ENV_PREFIX="$m2c_core"
HICREP_ENV_PREFIX="$m2c_hicrep"
DIFF_R_ENV_PREFIX=""
EOF
chmod 0644 "$m2c_site"
```

The profile contains software settings only. Project identity, inputs, biological design, outputs, scratch and resources belong in user configurations. It is trusted Bash and must not be writable by ordinary users. Do not add private credentials or dataset details to a publicly readable profile. Freeze this profile once users start analyses; use a new version for changed dependencies.

Create the release-specific command; keep its target pinned for reproducibility and resume:

```bash
m2c_launcher="$m2c_base/bin/microc2compare-0.1.0-shared-design-json-fix-20261008"
test ! -e "$m2c_launcher"
test ! -L "$m2c_launcher"
printf '#!/usr/bin/env bash\nset -Eeuo pipefail\nexport PYTHONDONTWRITEBYTECODE=1\nexec bash "%s/scripts/microc2compare.sh" "$@"\n' \
  "$m2c_src" > "$m2c_launcher"
chmod 0755 "$m2c_launcher"
test ! -e "$m2c_base/bin/microc2compare"
test ! -L "$m2c_base/bin/microc2compare"
ln -s "$(basename "$m2c_launcher")" "$m2c_base/bin/microc2compare"
"$m2c_launcher" --help
```

If a launcher/profile already exists, inspect it and use a new versioned name rather than overwrite a live deployment. The convenience name is for new projects; use the release-specific command/path for active runs and resume. No analysis command creates environments or uses sudo.

All accounts can use an absolute command immediately. Publish the convenience command through the site's existing `/opt/conda_envs/bin` mechanism if it is already on users' PATH. Otherwise the **administrator** can add one link in the usual system command directory. This requires no personal or system shell startup changes:

```bash
test ! -e /usr/local/bin/microc2compare
test ! -L /usr/local/bin/microc2compare
sudo ln -s "$m2c_base/bin/microc2compare" /usr/local/bin/microc2compare
```

Choose one publication location according to site policy; do not replace existing commands or create duplicate launchers. If a link already exists, inspect it and retain it until an explicitly reviewed promotion. Batch jobs and active runs should use the absolute versioned command. The convenience command runs this workflow only when called and does not activate Conda in other shells.

## 3. Start as an ordinary server user

After the maintainer activates the profile-bound 0.2.5 command, ordinary analysis requires your private project configuration:

```bash
microc2compare -c /path/to/your/project.conf --dry-run
microc2compare -c /path/to/your/project.conf
microc2compare -c /path/to/your/project.conf --analysis tads
```

Set `OUTDIR` and `SCRATCH_ROOT` to directories you own. No manual software-profile line is needed in the shell or project file. The optional workspace/cache setup below remains useful; its manual profile exports are needed for the separate upstream/demo/bootstrap helpers and older unbound commands, not for routine analysis through the new launcher.

**Each user, in a new login under their own account:** use the actual release/profile locations announced by the maintainer. The values below match this guide's example layout.

Run the shell command separately and wait for its prompt:

```bash
bash --noprofile --norc
```

Run this block from your own account. `$HOME` selects your home directory automatically:

```bash
set -Eeuo pipefail
unset HISTFILE BASH_ENV ENV PYTHONHOME PYTHONPATH R_HOME R_LIBS R_LIBS_USER
m2c_base="/opt/conda_envs/microc2compare"
m2c_src="/opt/conda_envs/microc2compare/releases/shared-design-json-fix-20261008/microc2compare-0.1.0"
m2c_site="/opt/conda_envs/microc2compare/site/shared-design-json-fix-20261008.conf"
source "$m2c_site"
umask 077
m2c_work="$HOME/microc2compare-work"
mkdir -p "$m2c_work/projects" "$m2c_work/validation" "$m2c_work/audit" \
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
/opt/conda_envs/microc2compare/bin/microc2compare-0.1.0-shared-design-json-fix-20261008 --help
```

Use a dedicated terminal for these session-only exports. Do not put them in `.bashrc` or run `conda init`/`conda config --set`; exit this terminal when finished. These variables select already installed software. This does not install packages, activate an environment globally or require a fresh copy of the source. A writable user temp directory is needed for configuration probes and Conda's command wrappers; see [Conda temporary-file guidance](https://docs.conda.io/projects/conda/en/latest/user-guide/configuration/temp-files.html). Record the actual Conda version and test it under the non-maintainer account; a historical version that tries to write inside its prefix is not accepted merely because it works for the owner. Fix that through a reviewed site installation, without granting ordinary users prefix-write access.

For a real project, export input drafts into a fresh directory in **your** workspace:

```bash
"$CONDA_EXE" run --no-capture-output -p "$CORE_ENV_PREFIX" \
  python "$m2c_src/scripts/prepare_microc2compare_inputs.py" \
  --samplesheet /ACTUAL/UPSTREAM/samplesheet.csv \
  --upstream-results /ACTUAL/UPSTREAM/results \
  --project my_project --reference-id mm39 \
  --write-drafts "$m2c_work/projects/my_project/config"
cp "$m2c_src/config/config_template.conf" \
  "$m2c_work/projects/my_project/config/project.conf"
printf '\nsource "%s"\n' "$m2c_site" >> \
  "$m2c_work/projects/my_project/config/project.conf"
```

Replace upstream paths/reference first and choose a new project destination if it exists. Read `READ_BEFORE_RUNNING.md` and [the handoff instructions](13-upstream-handoff.md). Edit the copied project configuration to match the actual draft filenames, reference registrations, biological tables, `PROJECT`/analysis/run IDs, relative `OUTDIR="../results"` and `SCRATCH_ROOT="../scratch"`. Source an appropriate resource preset **before** measured overrides. Keep optional pool/thinning/balancing disabled. Start with explicit descriptive comparisons; a profile with no differential R prefix does not support formal inference. Never edit files under the shared source tree for an analysis.

Run from any working directory using your private configuration:

```bash
m2c_project_conf="$m2c_work/projects/my_project/config/project.conf"
bash "$m2c_src/scripts/bootstrap_environments.sh" -c "$m2c_project_conf" --check-only
bash "$m2c_src/scripts/preflight_check.sh" -c "$m2c_project_conf"
/opt/conda_envs/microc2compare/bin/microc2compare-0.1.0-shared-design-json-fix-20261008 \
  -c "$m2c_project_conf" --dry-run
/opt/conda_envs/microc2compare/bin/microc2compare-0.1.0-shared-design-json-fix-20261008 -c "$m2c_project_conf"
```

## 4. Accept the installation from a second account

**An ordinary account, not the installation maintainer or root:** complete step 3's software/workspace setup, then check access. Run this before real projects and preserve the audit log. `namei -l` or `getfacl` can help the maintainer diagnose parent/ACL restrictions if installed.

```bash
test -r "$m2c_src/README.md"
test -r "$m2c_site"
test ! -w "$m2c_src"
test ! -w "$m2c_site"
test ! -w "$(dirname "$m2c_site")"
test ! -w /opt/conda_envs/microc2compare/bin
test ! -w "$CORE_ENV_PREFIX"
test ! -w "$HICREP_ENV_PREFIX"
test -x "$CONDA_EXE"
test -x "$PYTHON_BOOTSTRAP"
"$PYTHON_BOOTSTRAP" -c 'import sys; print(sys.version); assert sys.version_info >= (3, 11)'
"$CONDA_EXE" run --no-capture-output -p "$CORE_ENV_PREFIX" \
  python -c 'import cooler, cooltools, bioframe; print("Shared core imports OK")'
"$CONDA_EXE" run --no-capture-output -p "$HICREP_ENV_PREFIX" \
  python -c 'import hicrep; print("Shared HiCRep import OK")'
```

Expected: read/execute checks pass, ordinary users cannot write the shared software or replace the site profile/launcher, and both environment commands succeed. Root and the installation owner are unsuitable for these permission tests. A non-writable root can contain writable files, so inspect the complete source/prefix trees from this account without changing them:

```bash
find "$m2c_src" "$CORE_ENV_PREFIX" "$HICREP_ENV_PREFIX" -writable -print \
  > "$m2c_work/audit/unexpected-writable-software.txt"
test ! -s "$m2c_work/audit/unexpected-writable-software.txt"
```

Expected: the audit is empty and traversal succeeds. If scientific prefixes contain symlinks into an existing package cache, also review the targets' owner/ACL policy with the maintainer; this default `find` does not follow symlinks. Check the selected Conda/bootstrap installation and parent directories likewise. An intentionally privileged maintainer account is not an ordinary-user test. Do not diagnose permissions by creating files in an upstream prefix.

Now run a fresh six-sample smoke fixture in this account's workspace, entirely outside the release:

```bash
set -o pipefail
m2c_demo="$m2c_work/validation/shared-access-smoke-20261007"
"$CONDA_EXE" run --no-capture-output -p "$CORE_ENV_PREFIX" \
  python "$m2c_src/scripts/generate_demo.py" --out "$m2c_demo" \
  --core-prefix "$CORE_ENV_PREFIX" --raw-only
printf '\nsource "%s"\n' "$m2c_site" >> "$m2c_demo/project.conf"
sha256sum "$m2c_demo"/*.mcool > "$m2c_work/audit/shared-smoke-input-sha256.txt"
bash "$m2c_src/scripts/bootstrap_environments.sh" -c "$m2c_demo/project.conf" --check-only \
  2>&1 | tee "$m2c_work/audit/shared-smoke-software.log"
bash "$m2c_src/scripts/preflight_check.sh" -c "$m2c_demo/project.conf"
/opt/conda_envs/microc2compare/bin/microc2compare-0.1.0-shared-design-json-fix-20261008 \
  -c "$m2c_demo/project.conf" --dry-run
/opt/conda_envs/microc2compare/bin/microc2compare-0.1.0-shared-design-json-fix-20261008 \
  -c "$m2c_demo/project.conf" 2>&1 | tee "$m2c_work/audit/shared-smoke-run.log"
/opt/conda_envs/microc2compare/bin/microc2compare-0.1.0-shared-design-json-fix-20261008 \
  -c "$m2c_demo/project.conf" --resume
"$PYTHON_BOOTSTRAP" "$m2c_src/scripts/validate_outputs.py" \
  --config-json "$m2c_demo/results/demo/synthetic_demo/run001/00_metadata/config.json"
sha256sum --check "$m2c_work/audit/shared-smoke-input-sha256.txt"
```

Expected: six QC summaries, fifteen HiCRep pairs, a report, successful resume/output validation and unchanged input hashes. The fixture needs four CPUs/16,000 decimal MB by default; adjust both global and task budgets if your allocation is smaller. Repeat the full descriptive fixture from [installation step 10](15-server-installation-start.md#10-run-the-small-linux-smoke-test-then-the-full-descriptive-fixture) under this non-maintainer account, using its own new destination and omitting `--raw-only`. Both fixture modes must pass before accepting access to all enabled modules. Verify source manifest checksums again after the tests; no source changes should occur.

Test at least one additional ordinary account or representative access group for help/environment checks and a separate writable workspace. If users have different ACLs or mount access, validate each class; one second-account pass does not prove every account's policy. Confirm ordinary users can access the **software**, while unauthorized access to private project data remains denied.

## 5. Coordinate resources across users

`TOTAL_CPUS`, `MAX_MEMORY_MB` and concurrency caps apply to **one run**. Preflight observes current availability but cannot atomically reserve RAM/CPUs across independent users. The reported 500 GB/70 physical cores do not become separate budgets for each account. Shared installation does not enforce a server-wide quota, and run locks prevent output collisions rather than aggregate overcommit.

Use the existing scheduler where available: every job requests the CPUs/memory matching its project budget and runs inside that allocation, with site scratch paths selected normally. The shell workflow needs no Snakemake/Nextflow to run inside a scheduler job. Batch jobs should use absolute release/profile paths because login PATH initialization may not run. On a server without a scheduler, agree a central submission/resource policy with the administrator and coordinate heavy runs before launch; until one exists, run one heavy real-data pilot at a time. Do not advertise the 60-CPU/380-GB production candidate as a simultaneous per-user default.

Begin with the demo budget, then choose a measured per-job pilot allocation. Account for other workflows, disk/inode quotas, I/O load and node-local scratch. See [resource tuning](07-server-resources.md). Software availability for all users does not imply unlimited concurrent access to the hardware.

## 6. Updates, rollback and troubleshooting

Install new source/releases, clone or fallback prefixes, profiles and launchers alongside old versions. Test as maintainer and ordinary user, retain the old paths, and switch the convenience command only after acceptance. Keep active runs pinned to their original release/profile/prefixes through resume. Do not upgrade a reused upstream prefix or edit a versioned profile in place. A new published source/environment combination still needs its applicable numerical, Linux and statistical validation gates.

| Problem | Action |
|---|---|
| Works only for installer | Check all parent traversal permissions, ACLs, Conda/bootstrap paths and activation hooks; repeat the non-maintainer test |
| `microc2compare` not found | Use the absolute versioned command; confirm the site module/login PATH policy or add shared `bin` in this shell |
| Demo/config cannot be written | Use your private workspace, not the shared checkout |
| Conda wrapper/cache needs a writable path | Check user `TMPDIR`/cache and actual Conda behavior; never grant prefix-write access as a workaround |
| Upstream dataset unreadable | Obtain that dataset's authorized read access from its owner; software installation cannot grant data access |
| Another user has a run lock | Select your own results/analysis/run path; do not delete another account's lock |
| Several users exhaust RAM/I/O | Reduce per-job budgets and use the site's scheduling policy; per-run preflight is not global admission |
| Shared profile changed during an active run | Restore the original versioned software settings or start a new reviewed run; do not edit completion fingerprints |

This document describes the required deployment and acceptance steps. No administrator commands, shared server permissions or second-account Linux tests were executed during local preparation.
