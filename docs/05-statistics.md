# Statistical interpretation and gates

The report's [differential loop volcano plots](22-loop-volcano-plots.md) show model log2 allocation effects against biological P values, with colours determined by full-family BH FDR. They never substitute caller enrichment statistics for biological evidence. Descriptive comparisons have no volcano; preliminary/development limitations remain labelled. Optional magnitude and label settings change display only.

One manifest row is one biological unit/condition observation. Technical libraries, lanes, pools and thinning seeds are never independent replicates. Distinct IDs alone do not prove independence. Register actual experimental units and deduplication history before interpreting any contrast.

## Eligible comparisons

| Requested mode | Requirement and output |
|---|---|
| `descriptive` | At least one observation per condition; effects only, biological P/FDR are NA |
| `preliminary` | Explicit request, at least two independent units per condition, full rank, residual df ≥2, independent loop universe and matching calibration |
| `formal` | At least three units per condition, full rank, residual df ≥2, independent loop universe and matching calibration |

`~condition` with two unpaired samples per condition has two residual df and can be preliminary. `~subject_id+condition` with two matched pairs has only one residual df and fails this policy. Three pairs leave two residual df. Completely confounded assay/batch/condition terms fail rank checking. Cross-assay requests are descriptive assay-sensitivity comparisons. Failure never silently downgrades inference.

Formula terms are additive only. Subject/batch/assay/condition labels are factors, including numeric-looking labels. Other numeric metadata becomes a continuous covariate and text becomes a factor, with the control condition as the reference. The contrast is case minus reference. Record the exact exported design matrix/contrast and inspect leverage, confounding, outliers and metadata accuracy in the biological pilot. Additional covariate designs need their own calibration domains; the shipped simulation runner supports unpaired condition-only and complete-pair subject-plus-condition designs with explicitly selected unit counts.

## Formal loops: relative allocation

Freeze an independently discovered external BEDPE universe with inspectable independence evidence before inspecting the comparison outcomes. Same-data dots/consensus and condition-blind pools do not satisfy this requirement.

For feature `i` in biological sample `s`, measure the single raw 10-kb anchor-pair pixel `Y_is`. Define one scalar `L_s` as all unique raw eligible same-view cis contacts in the common mask and inclusive 50-kb–2-Mb range across all registered views. Fit:

```text
log(mu_is) = log(L_s) + X_s beta_i
```

The R adapter uses edgeR robust dispersion estimation and quasi-likelihood fitting/testing with library sizes `L_s` and normalization factors exactly one. It uses no TMM, balancing weight, O/E, subtraction, candidate-count column sum or pseudoreplicate. This estimates **relative contact allocation**, not absolute contacts or focal enrichment. A global redistribution can change this estimand even if a candidate's absolute contacts do not change.

Prespecified filtering retains complete eligible features with total raw count ≥20 and `mean(Y_is/L_s)*min(L_s)≥5`; at least 1,000 retained independent features are initially required for dispersion learning. These are pilot policies requiring calibration, not universal guarantees. BH spans all retained tested genome-wide features in the declared comparison family, never individual chromosomes or selected significant subsets.

Outputs retain raw counts, exposure/design inputs, log2 allocation effect, abundance, QL statistic, dispersion, degrees of freedom, P value and BH FDR. The approximate QL-information Wald standard error and pointwise interval are clearly labeled approximate; the QL F-test supplies the P value. These intervals do not claim simultaneous coverage and require numerical/coverage validation with the selected edgeR version. Filtering/exclusion reasons remain explicit.

## Formal insulation: fixed grid

Use finite comparable log2 scores on the common fixed 10-kb grid at the prespecified primary 100-kb window. Eligibility/mask/coverage rules precede condition testing. The adapter directly uses limma `lmFit`, `contrasts.fit`, `eBayes(robust=TRUE, trend=TRUE)`; it does not use voom on continuous insulation scores.

BH spans all tested genome-wide eligible bins in this primary family; boundaries are annotations afterward. Export moderated SE, t, df, P, BH FDR and pointwise 95% intervals. Additional insulation windows are separate descriptive sensitivity families in v0.1. Merged significant intervals are descriptive and have no inherited regional FDR. Missing/non-minimum prominence remains NA.

Insulation is relative to its algorithmic within-view normalization. Spatial correlation, small sample sizes, depth effects and heavy-tailed/outlier scores require calibration and residual inspection. Biological effects and coverage artifacts must be separated in the pilot.

## Calibration and development labels

Formal/preliminary runs require a passing artifact matching the exact source, environments, module, design terms, biological group sizes, subject/observation counts, design rank/df and filter/model recipe. A three-replicate calibration cannot certify a two-replicate preliminary design. Empty, mismatched, failed or underpowered evidence blocks preflight. Each primary scenario requires ≥200 independent simulated datasets and upper 95% Monte Carlo mean-FDP bound ≤0.075 at nominal FDR 0.05. [The protocol](../validation/calibration-protocol.md) specifies simulation truth, uncertainty and pairing. Select the actual group size with `--units-per-condition`; two matched pairs fail the df policy even in calibration.

`ALLOW_UNCALIBRATED_DEVELOPMENT=true` is an explicit method/fixture development exception. Every comparison table and run report is labeled `validation_status=uncalibrated`; it cannot satisfy the release gate. No passing artifact is shipped. A passing simulation key supports only its declared design domain, not arbitrary covariates or real-data validity. Biological negative controls and public-data validation remain separate mandatory review gates.

```bash
"$CONDA_EXE" run --no-capture-output -p "$CORE_ENV_PREFIX" \
  python scripts/calibrate.py --config-json /path/to/run/00_metadata/config.json \
  --module loops --design unpaired --datasets 200 \
  --out "$HOME/microc2compare-work/calibration/loops"
# Repeat for insulation and the separate paired design domain if it will be used.
# Register reviewed existing evidence with scripts/register_calibration.py;
# see docs/19-r-environment-and-inference.md for complete commands.
```

The simulation procedure exercises the actual shipped R adapters. It does not create independent biological observations or replace negative controls. The runner records power and FDP; further effect-bias and interval-coverage review is required before endorsing approximate uncertainty fields.

Primary method references: [edgeR](https://bioconductor.org/packages/edgeR/), [limma](https://bioconductor.org/packages/limma/). Caller significance, biological FDR, aggregate enrichment and descriptive effect magnitudes must not be conflated.
