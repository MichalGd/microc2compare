# Differential loop volcano plots and amended candidate installation

For a new installation from repository source, use [the GitHub installation guide](26-install-from-github.md). ZIP transfer sections here document earlier private deliveries; those artifacts are not included in Git.

**Latest acceptance correction:** 0.2.3 includes these features and corrects raw eigenvector text comparisons. For the reported `compartments/C1.tsv differ` failure, follow [numerical agreement and resume](23-compartment-agreement.md). The plot documentation/settings below remain current; the 0.2.2 installation steps retain release history.

Version **0.2.2** includes the TAD/module additions, the dependency-free bootstrap correction and differential loop volcano plots. Install this candidate directly; installing 0.2.1 first is unnecessary. Existing scientific environments and imported matrices need no changes. The global command stays on its previously accepted release until candidate acceptance passes.

## What the plot means

Each point is a retained loop tested by the biological-replicate model. The horizontal axis is edgeR `logFC`: log2 relative contact allocation in the case condition divided by the reference condition. Positive values mean more allocation in the case; negative values mean less. An effect of 1 corresponds to twice the relative allocation. This does not measure an absolute genome-wide change in contact abundance.

The vertical axis is `−log10(PValue)` from the differential model: higher points have smaller biological P values. Colours use the model's Benjamini–Hochberg FDR across its complete retained testing family. FDR is the estimated proportion of false discoveries among results selected at that threshold. Loop-caller enrichment P/q values are never used. There is no horizontal line presented as an FDR cutoff: the raw P-value cutoff corresponding to an adjusted FDR depends on the complete testing family.

The plot labels the direction and inference status. Preliminary results remain preliminary. Explicitly permitted development analyses are marked **UNCALIBRATED DEVELOPMENT** and cannot support formal biological conclusions. A calibrated label reflects the existing calibration gate; the plot itself provides no calibration or biological validation.

The current descriptive synthetic installation demo has no biological P values/FDR and generates **no volcano plot**. Its HTML/status JSON explains this. Inferential plots require the checked R/edgeR environment, independent frozen loop candidates, sufficient biological replicates, an appropriate design and the existing validation gates. See [R setup](19-r-environment-and-inference.md) and [statistics](05-statistics.md). Synthetic model-result fixtures in local plot tests check rendering, not actual R execution/calibration.

## Inputs, outputs and settings

The report reads `06_comparisons/ID/loops.tsv` and `parameters.json` only when both are declared outputs of the currently valid differential completion. It uses configured loop comparisons only. No input matrix is copied or balanced and no statistical model is rerun by the plotting helper. Stale/uncommitted files cannot authorize a plot.

Required inferential columns are `feature_id`, `logFC`, `PValue`, `FDR`, `retained`, `inference_status` and `validation_status`. Only retained rows with finite effects and P values/FDR between 0 and 1 are plotted. Filtered rows and invalid statistics are counted. A zero P value caused by numeric underflow uses the smallest positive normal float for display only; the plotted-points TSV preserves zero and records `pvalue_display_clipped=true`. Original model results remain untouched.

| Output under the run directory | Meaning |
|---|---|
| `08_reports/index.html` | Embedded plots, direction/status and explanations for unavailable plots |
| `08_reports/figures/volcano.ID.png` | Raster image for viewing/slides |
| `08_reports/figures/volcano.ID.svg` | Scalable axes/text; points rasterized to keep large plots manageable |
| `08_reports/volcano/ID.points.tsv` | All usable plotted rows, original effects/P/FDR, height, FDR flag, colour group and zero-P display flag |
| `08_reports/volcano/status.json` | Reasons, input/plotted/filtered/invalid counts, direction, inference/validation and thresholds |

PNG/SVG/point exports exist only for available plots. The status JSON remains available for enabled comparisons even when no plot is appropriate. Empty/invalid-only results have an explicit reason. The HTML links current exports only; earlier files are not presented as current products.

Project `.conf` settings:

```bash
PLOT_LOOP_VOLCANO="true"
VOLCANO_MIN_ABS_LOG2_EFFECT="0"
VOLCANO_LABEL_TOP_N="5"
```

Set `PLOT_LOOP_VOLCANO=false` to disable plots. `VOLCANO_MIN_ABS_LOG2_EFFECT` is a nonnegative magnitude for colouring/labels only: all usable points remain plotted and no P value, FDR or testing family changes. `VOLCANO_LABEL_TOP_N` labels up to that many highlighted loops with the smallest model P values; 0 disables labels, maximum 20. Significance uses the saved differential `FDR_THRESHOLD` (default 0.05). Changing that scientific threshold affects comparison signatures; the three plot settings affect report signatures only.

## Launch separately or with the pipeline

Set your actual project configuration and source path in server Bash:

```bash
m2c_src="/opt/conda_envs/microc2compare/releases/volcano-20261008/microc2compare-0.2.2"
m2c_conf="$HOME/microc2compare-work/projects/experiment/project.conf"
bash "$m2c_src/scripts/microc2compare.sh" -c "$m2c_conf" --analysis loop_comparison --dry-run
bash "$m2c_src/scripts/microc2compare.sh" -c "$m2c_conf" --analysis loop_comparison --resume
```

Replace the illustrative `experiment` path with your real config. `loop_comparison` enables required quantification and the report, preserving statistical gates. Omit `--analysis` for the configured whole pipeline. After display changes, resume **the same original selection** without stage bounds: preflight refreshes the configuration, compatible scientific completions are reused and the report regenerates. Plot settings do not force model recomputation.

After a successful run/resume with the current settings, you may explicitly rerender only the report:

```bash
bash "$m2c_src/scripts/microc2compare.sh" -c "$m2c_conf" \
  --analysis loop_comparison --from-stage report --until-stage report --rerun-stage report --resume
```

This requires current prerequisite completions, including preflight. If settings, source, inputs, model parameters or selection changed, resume the original selection without stage bounds first. `--analysis report` enables configured prerequisites and cannot bypass them. Preserve historical parameter comparisons with [a variant/new run ID](17-analysis-modules.md#reports-variants-and-errors).

## Transfer: Windows PowerShell

Run on Windows at a prompt containing `PS`, not on the server:

```powershell
$m2c_delivery = "D:\path\to\microc2compare\.local"
scp -P 22 "$m2c_delivery\microc2compare-0.2.2-volcano-8oct26.zip" "$m2c_delivery\microc2compare-0.2.2-volcano-8oct26.zip.sha256" "$m2c_delivery\install-microc2compare-0.2.2.py" "USERNAME@SERVER:/opt/conda_envs/microc2compare/incoming/"
```

## Verify and install: server Bash

Use the dedicated Bash/private cache settings from [candidate installation](20-candidate-upgrade.md#2-start-a-dedicated-linux-session). Assign variables before using them:

```bash
m2c_base="/opt/conda_envs/microc2compare"
m2c_work="$HOME/microc2compare-work"
m2c_bootstrap="/opt/miniconda/bin/python"
cd "$m2c_base/incoming"
sha256sum --check microc2compare-0.2.2-volcano-8oct26.zip.sha256
```

Continue only if the checksum reports `OK`:

```bash
read -r m2c_sha m2c_archive_name < "$m2c_base/incoming/microc2compare-0.2.2-volcano-8oct26.zip.sha256"
"$m2c_bootstrap" "$m2c_base/incoming/install-microc2compare-0.2.2.py" \
  --archive "$m2c_base/incoming/microc2compare-0.2.2-volcano-8oct26.zip" \
  --sha256 "$m2c_sha" --base "$m2c_base" --release-id volcano-20261008
m2c_src="$m2c_base/releases/volcano-20261008/microc2compare-0.2.2"
bash "$m2c_src/scripts/microc2compare.sh" --help
```

The installer verifies the manifest, refuses an existing release ID and leaves earlier source/envs/launchers untouched. If this ID already exists, inspect the earlier installation result; do not overwrite it.

## Resume the failed acceptance fixture: server Bash

The original maps/configurations/hash list are retained. The helper verifies their hashes before scientific execution and does not regenerate inputs or reset identities:

```bash
m2c_site="$m2c_base/site/shared-design-json-fix-20261008.conf"
m2c_accept="$m2c_work/validation/modules-candidate-20261008"
set -o pipefail
bash "$m2c_src/scripts/validate_module_demo.sh" \
  --site "$m2c_site" --out "$m2c_accept" --resume \
  2>&1 | tee "$m2c_work/audit/modules-candidate-volcano-20261008.log"
m2c_accept_status=${PIPESTATUS[0]}
printf 'Candidate acceptance exit code: %s\n' "$m2c_accept_status"
```

Success requires `Full/selected/variant acceptance passed`, exit code zero and unchanged original hashes. Source/schema changes can correctly recompute earlier tasks. This descriptive demo checks the explanation for absent plots; local explicit model-result fixtures check rendering. Actual inferential plot acceptance follows R numerical/design calibration acceptance. Retain logs if any step stops.

Repeat ordinary-account execution in that account's chosen private workspace. Only software/source are shared; configurations, results, scratch and caches remain user-owned. GitHub publication remains deferred.

## Publish only after candidate acceptance

After acceptance, ordinary-account checks and report review, the maintainer creates the versioned launcher. This leaves the existing shared alias in place:

```bash
"$m2c_bootstrap" "$m2c_src/scripts/publish_shared_command.py" --source "$m2c_src" --base "$m2c_base"
"$m2c_base/bin/microc2compare-0.2.2-volcano-20261008" --help
```

Then activate the accepted source explicitly:

```bash
"$m2c_bootstrap" "$m2c_src/scripts/publish_shared_command.py" --source "$m2c_src" --base "$m2c_base" --activate
hash -r
microc2compare --help
```

The existing `/usr/local/bin/microc2compare` symlink needs no changes. Keep earlier versioned launchers/releases. The [earlier publication explanation](20-candidate-upgrade.md#6-publish-the-accepted-command) describes the same process; the commands here use the amended version.

## Troubleshooting

- **No volcano in the demo:** expected for `inference=descriptive`; inspect `volcano/status.json`.
- **Missing/stale differential output:** resume the original selection; do not hand edit completion records.
- **All points unavailable:** inspect `retained`, invalid-statistic counts and filtering evidence in the original table.
- **Missing R/calibration:** follow the R/design guide; caller P values cannot substitute for biological statistics.
- **Downloaded HTML unchanged:** download the complete updated run folder, including figures/volcano data; reopen the local HTML.
