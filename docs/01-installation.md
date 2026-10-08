# Linux installation

For the complete numbered path from workstation transfer through existing-environment checks, isolated fallback creation, Linux fixtures and upstream input preparation, start with [step-by-step server installation](15-server-installation-start.md). This page records the environment/lock policy in more detail.

The required server deployment is a maintainer-controlled shared release accessible to all server users, with separately writable user project/results/scratch/cache directories. Follow [shared server installation](16-shared-server-installation.md) for provisioning, the versioned site profile/command and non-maintainer acceptance. Do not install shared software in a private home or generate demos/configurations inside the published checkout. Reuse requires both scientific compatibility and approved read/execute access for intended accounts; never broaden upstream/home permissions to force reuse.

The target is a Linux x86-64 server. Bash 4.4+ supplies `mapfile`, arrays and process-group handling; bootstrap Python 3.11+ parses manifests/configuration without scientific dependencies. Local verification used Bash 5.3.9 and Python 3.11.9 on Windows. Linux runtime support must pass the fixture/CI before release.

## Place the checkout and project files

Transfer the reviewed public archive or clone the later published repository into an owned directory. Do not install as root merely to use the workflow. Keep environments, immutable references and projects in separate owned locations. The following uses example writable `/opt/conda_envs/microc2compare` directories; change them to approved local paths:

```bash
m2c_src="/opt/conda_envs/microc2compare/releases/shared-design-json-fix-20261008/microc2compare-0.1.0"
m2c_project="$HOME/microc2compare-work/projects/pilot"
umask 077
mkdir -p "$m2c_project/config"
test ! -e "$m2c_project/config/project.conf" && \
  cp "$m2c_src/config/config_template.conf" "$m2c_project/config/project.conf"
```

Edit the copied file. Set stable identifiers and absolute manifest/reference paths. Source the checkout's `config/resources_pilot.conf` before any project-specific resource overrides. Set `OUTDIR` and `SCRATCH_ROOT` to owned paths with measured free space. Do not source an upstream project configuration; register its relevant provenance explicitly.

## Check selected immutable prefixes

Before real manifests are complete, use the guide's software-only `setup.conf` with `RUN_DIFFERENTIAL=false`. With differential analysis enabled, environment discovery reads comparisons/biological metadata; copying the project template with nonexistent sample tables is therefore not a valid first bootstrap configuration.

The recommended deployment clones the existing Linux upstream core once into a new shared prefix, using a central manager outside homes. Follow installation step 6 for inventory, guarded cloning, target-only additions and checks; step 7 is the clean fallback. Direct reuse remains technically supported when explicitly chosen and compatible, but a separate clone keeps future package changes independent of upstream work. Discover actual paths from the server administrator's installation record. There is no assumption that `/opt/conda/envs/microc2tracks` exists. The commands below inspect only:

```bash
command -v conda
conda list -p /actual/upstream/prefix --json
conda run --no-capture-output -p /actual/upstream/prefix python --version
```

Set `ENVIRONMENT_MODE="reuse"`, `CONDA_EXE`, `CORE_ENV_PREFIX`, `HICREP_ENV_PREFIX`, and, for inferential requests, `DIFF_R_ENV_PREFIX`. HiCRep can use the same checked prefix if already compatible, or an isolated project prefix. Differential R is not required by descriptive comparisons.

```bash
bash "$m2c_src/scripts/bootstrap_environments.sh" -c "$m2c_project/config/project.conf" --check-only
```

This checks required imports/APIs and small raw/weighted/HiCRep numerical contracts. It does not certify real-data validity or update packages. Preflight records package inventories, Python and R runtime fingerprints. `PYTHONNOUSERSITE`, disabled user R startup files and nested-thread caps prevent accidental user-package/thread contamination.

## Create missing owned prefixes explicitly

If a reused prefix fails checks, point the failed environment at a **new** owned location. Existing non-Conda directories are rejected; existing Conda prefixes are checked without mutation.

```bash
# After editing prefixes in the project configuration:
# Maintainer only: provision new shared prefixes, not one copy per user.
bash "$m2c_src/scripts/bootstrap_environments.sh" -c "$m2c_project/config/project.conf" --create-missing
bash "$m2c_src/scripts/bootstrap_environments.sh" -c "$m2c_project/config/project.conf" --check-only
```

Core candidates are Cooler 0.10.4, cooltools 0.7.1, bioframe 0.7.2 and Python 3.11; `environment.yml` pins the principal scientific packages. HiCRep 0.2.6 uses `envs/hicrep.yml` and a separate pip installation performed by Conda environment creation. The R definition uses R 4.5, edgeR 4.8.2 and limma 3.66.0 after the original R 4.6/edgeR 4.10.5/limma 3.68.5 combination failed to resolve. The revised solve selected R 4.5.3. Solve evidence is in [environment verification](../validation/environment-solves.json).

Principal version pins and successful cross-platform solves are not complete tested locks. The HiCRep dry solve verifies Conda dependencies; its pip package was installed/tested locally on Windows, not by the Linux dry solve. Installation still resolves transitive dependencies unless a lock is explicitly configured. No solve or install occurs during analysis. Never run `conda update` on an upstream shared prefix to satisfy this workflow.

## Exact locks after Linux validation

`ENVIRONMENT_MODE="locked"` requires explicit lock paths and matching `LOCK.inventory.json` package digests. Bootstrap uses `conda create -p PREFIX --file LOCK` directly; missing locks cannot trigger a fallback solve. The initial release contains no fabricated tested locks. Follow [lock instructions](../envs/locks/README.md) after the full Linux fixture passes. HiCRep additionally requires `LOCK.pip.txt` with reviewed hashed wheel requirements, installed with `--require-hashes --no-deps` in a newly created owned prefix. Missing pip records fail instead of a silent solve. Until this is prepared and verified, use checked immutable prefixes in `reuse` mode.

Record OS/glibc, architecture, Conda solver, complete package JSON, executable paths and checksums. Validate a recreated prefix before registering its inventory. A lock inventory is a drift check, not evidence of calibration. Exported absolute prefix fields must be canonicalized consistently by the supplied fingerprint implementation; compare inventories on the same Conda command format.

## Validate before the server pilot

Generate the synthetic fixture in an empty directory, set checked core/HiCRep prefixes, and use the demo resource preset. The four-CPU/16,000-MB demo budget must be lowered on a smaller allocation; override every enabled task memory reservation accordingly. Execute the full descriptive fixture, validate completion hashes, resume it, and confirm input SHA256 values are unchanged. Then inventory the real server and proceed with [the measured pilot](07-server-resources.md).

No remote installation, service changes or production datasets were executed during local development.
