# R environment inspection and statistical acceptance

**Current server status:** the user reported passing 0.2.4 R numerical checks, non-owner R package access and shared activation using the existing `/opt/conda_envs/rna_ends2tracks-0.1.0a12.post1` prefix. No additional R installation is required for that accepted combination. Statistical calibration remains pending. A [profile-bound 0.2.5 launcher](25-automatic-software-profile.md) can select the same checked prefix automatically; separate administrative helpers retain explicit settings.

R is optional for descriptive analysis, TAD events and reports. It is required for preliminary/formal loop allocation and fixed-grid insulation comparisons. Finding R or DESeq2 alone is insufficient: the shipped adapters use edgeR, limma, statmod and jsonlite. Installation and numerical acceptance do not establish biological/statistical calibration.

If an older 0.2.3 loop numerical check failed with `vapply ... result is length 0`, install [the 0.2.4 R adapter correction](24-r-loop-compatibility.md) before continuing. Local actual-R tests ran with edgeR 4.4.2; the user subsequently reported acceptance with edgeR 4.4.0 on this server.

## 1. Inspect existing shared environments first

These are **Linux Bash** commands after the dedicated-session setup. They list existing prefixes without Conda environment registry cleanup, installation or package changes. A path being listed does not prove compatibility.

```bash
m2c_src="/opt/conda_envs/microc2compare/releases/r-ql-fix-20261008/microc2compare-0.2.4"
m2c_bootstrap="/opt/miniconda/bin/python"
"$m2c_bootstrap" "$m2c_src/scripts/find_r_environments.py"
```

Choose a candidate shared prefix outside users' homes. It must remain immutable while analyses use it. Inspect its ownership/read permissions and coordinate with its maintainer. Reusing an environment never authorizes modifying an ATAC/RNA environment to add missing packages.

## 2. Run read-only package and numerical checks

Use your checked core Python to run the acceptance helper. It runs the actual shipped R adapters on small synthetic tables and saves evidence in a **new private directory**, leaving the environment unchanged. Replace the candidate path first:

```bash
m2c_site="/opt/conda_envs/microc2compare/site/shared-design-json-fix-20261008.conf"
source "$m2c_site"
m2c_r_candidate="/REPLACE/WITH/EXISTING/SHARED/R/PREFIX"
m2c_r_audit="$HOME/microc2compare-work/validation/r-candidate-check-20261008"
"$CONDA_EXE" run --no-capture-output -p "$CORE_ENV_PREFIX" \
  python "$m2c_src/scripts/validate_r_environment.py" \
  --conda "$CONDA_EXE" --prefix "$m2c_r_candidate" --out "$m2c_r_audit"
```

Expected: package/session evidence, raw model inputs, loop/insulation result TSVs, per-command logs and `acceptance.json` with `status=passed_numerical_checks_only`, `calibration_status=not_run`. Checks cover effect direction, null effects, exposure scaling, missing-score exclusion, uncertainty fields and whole-tested-family BH adjustment. A failed check exits nonzero; inspect its log and retained failed status. The helper refuses an existing audit directory.

If no suitable prefix passes, the maintainer can create a new versioned owned prefix from `envs/diff-r.yml`, using the existing shared Conda manager and the dedicated private package cache. No user's ordinary analysis command performs this installation:

```bash
m2c_r_candidate="/opt/conda_envs/microc2compare/envs/diff-r-0.2.4-20261008"
test ! -e "$m2c_r_candidate"
umask 022
"$CONDA_EXE" env create -p "$m2c_r_candidate" -f "$m2c_src/envs/diff-r.yml"
umask 077
```

Execute the creation only if the preceding existence check succeeds. Restore `umask 077` afterward, then repeat numerical acceptance in a new audit directory. Run in a dedicated shell with the installation guide's `CONDA_REGISTER_ENVS=false`, private cache/temp and no shell startup edits. Keep the accepted core and HiCRep prefixes unchanged. Confirm ordinary-account read/execute access and reject group/world-writable published software.

Record an exact package inventory/lock and test reconstruction before endorsing an immutable prefix. Put its path in a **new versioned software profile** using `DIFF_R_ENV_PREFIX`; retain the earlier descriptive profile. Software profiles contain only software paths/mode/locks, not shared output directories or datasets. Source the chosen profile at the end of each private project configuration.

## 3. Prepare a valid comparison and calibration configuration

`comparisons.tsv` declares `descriptive`, `preliminary` or `formal` per row. `RUN_DIFFERENTIAL=true` only enables the adapters; it does not choose formal inference. Formal comparisons require at least three genuine independent units per condition and full-rank designs with at least two residual degrees of freedom. Preliminary requests require at least two independent units and their own passing design/calibration. Matched subjects require the declared paired design; technical libraries/lanes are not units.

Formal loops require an independently discovered frozen feature registry with evidence, matched assembly/grid, and enough retained candidates for the prespecified model (initial floor 1,000). Same-data consensus or a condition-blind pool does not establish independence. Formal insulation uses the common fixed primary grid. See [statistics](05-statistics.md) for estimands, filtering, paired designs and exclusions.

To prepare calibration inputs, create a **separate development configuration/run ID** with the intended inferential design, checked R profile and `ALLOW_UNCALIBRATED_DEVELOPMENT=true`. This is a labeled method-development exception, never a production result. Run the intended module with that explicit configuration to produce its `00_metadata/config.json`; do not weaken replicate/feature requirements or silently relabel an existing result. Use that snapshot below and a user-owned registry path.

## 4. Run and review the prespecified simulation suites

These potentially substantial computations are separate from installation. Reserve a suitable per-run CPU/memory budget and use a private output directory. The example is for an unpaired three-unit design; use the actual intended group size and repeat the supported paired design if needed.

```bash
m2c_cal_config="$HOME/microc2compare-work/projects/calibration/results/demo/analysis/dev-run/00_metadata/config.json"
m2c_cal_out="$HOME/microc2compare-work/validation/calibration-loops-20261008"
"$CONDA_EXE" run --no-capture-output -p "$CORE_ENV_PREFIX" \
  python "$m2c_src/scripts/calibrate.py" --config-json "$m2c_cal_config" \
  --module loops --design unpaired --units-per-condition 3 --datasets 200 \
  --out "$m2c_cal_out"
```

Replace the configuration path with the actual development run snapshot. Repeat for insulation and each supported intended design domain. The runner evaluates multiple primary scenarios with at least 200 independent datasets each. Passing evidence must match source, environments, design/group sizes and model/filter recipe. Review false discoveries, power, effect bias, interval coverage and real-data/negative-control behavior; simulation acceptance alone does not certify arbitrary biological data.

Register only reviewed passing evidence. The registration helper reuses the actual preflight gate, checks the current source/environments/design/recipe, and avoids rerunning the simulations:

```bash
m2c_cal_registry="$HOME/microc2compare-work/validation/calibration_registry.json"
"$CONDA_EXE" run --no-capture-output -p "$CORE_ENV_PREFIX" \
  python "$m2c_src/scripts/register_calibration.py" --config-json "$m2c_cal_config" \
  --results "$m2c_cal_out/calibration_results.json" --registry "$m2c_cal_registry" \
  --comparison-id treated_vs_control --module loops
```

Replace the comparison ID with the actual intended contrast. Repeat with insulation evidence as needed. Configure `CALIBRATION_REGISTRY` to the owned registry, and use a **new production run ID** with `ALLOW_UNCALIBRATED_DEVELOPMENT=false`. Failed/mismatched evidence cannot be registered. Do not hand-edit registry keys/status or reuse evidence from another source/recipe/design.

## 5. Execute and interpret

```bash
m2c_inference_conf="$HOME/microc2compare-work/projects/pilot/config/formal.conf"
bash "$m2c_src/scripts/microc2compare.sh" -c "$m2c_inference_conf" --analysis comparisons --dry-run
bash "$m2c_src/scripts/microc2compare.sh" -c "$m2c_inference_conf" --analysis comparisons
```

Preflight must pass every design/input/calibration requirement. Results live in `06_comparisons/ID/`. Inspect effects, exclusions, diagnostics, uncertainty and `inference_status`/`validation_status`, not just FDR. Loop tests estimate relative contact allocation; insulation tests concern continuous relative scores. TAD events, domain geometry, compartments and aggregate plots remain descriptive and do not acquire event/domain FDR from nearby significant bins.

The correction has been numerically tested in an isolated local Windows R runtime; those versions and scope are documented in [the correction guide](24-r-loop-compatibility.md). No passing Linux R acceptance or calibration artifact is shipped with this candidate. Record what actually ran in your private audit before enabling inference.
