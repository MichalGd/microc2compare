# Analysis catalogue: inputs, commands, settings and outputs

**0.3.0 additions:** `chromosight` and `annotation` are independently selectable. Loop discovery/quantification/comparison and APA support complete resolution lists, while ADA uses `ADA_RESOLUTION` and saddle uses `COMPARTMENT_RESOLUTION`. Autosome selection is the default. For new commands, inputs, settings, per-resolution outputs and interpretation, use [the fine-resolution module table](28-fine-resolution-loops.md#5-separate-modules-and-outputs) and [anchor annotation](29-loop-anchor-annotation.md). Primary-resolution aliases preserve older output paths below. The published command handles software settings and private caches internally; no separate user `source` or clean-shell session is required.

Use this guide to choose an analysis and find its results. A contact matrix records how often two genomic bins interact. A bin is a fixed-length genome interval. A biological sample is one independent animal, donor or culture, with its technical libraries consolidated upstream.

Every analysis reads original contact counts and, where required, existing upstream balancing weights. No command below installs packages or creates a copy of an imported matrix. Results are private/user-managed: `OUTDIR/PROJECT/ANALYSIS_SET_ID/RUN_ID`; temporary work uses `SCRATCH_ROOT`. Both paths are configurable. Relative paths in `project.conf` resolve against that file's directory.

## Start, select and resume

These are Linux Bash commands in your normal login. The accepted shared launcher supplies software and runtime settings automatically. Set only your project-specific configuration path:

```bash
m2c_conf="$HOME/my_microc_project/config/project.conf"
microc2compare --list-analyses
microc2compare -c "$m2c_conf" --analysis tads --dry-run
microc2compare -c "$m2c_conf" --analysis tads
microc2compare -c "$m2c_conf" --analysis tads --resume
microc2compare -c "$m2c_conf" --validate
```

The first execution needs a new run directory. An existing run requires `--resume`. Keep the same selection when resuming. If you deliberately change selection, registration/dependent results may need refreshing; the runner checks compatibility. Necessary prerequisites are enabled automatically and run in their usual order. Disabled analyses create no placeholder scientific files. A report is produced unless `RUN_FINAL_REPORT=false`.

Multiple names are comma-separated, for example `--analysis compartments,saddle`. Repeating `--analysis` combines selections. Omit `--analysis` for the configured full pipeline. `--analysis report` retains the configured science and refreshes/reuses its prerequisites; it is not a request to present unchecked stale results. Existing `--from-stage`, `--until-stage` and `--rerun-stage` controls remain available. Starting in the middle still requires current prerequisites; start from preflight after changing settings.

## Common inputs

All scientific modules need the biological samples table, technical membership table, registered chromosome sizes/views, and compatible matrix grids. See [exact input contracts and example records](02-inputs.md) and [upstream preparation](13-upstream-handoff.md).

| File | Expected content |
|---|---|
| `samples.csv` | `sample_id,project,condition,assay,biological_replicate,subject_id,reference_id,matrix_path,provenance_path`; paths identify consolidated biological maps and processing evidence |
| `technical_membership.tsv` | Sample/library/source relationships, source kind/checksum, deduplication group/scope; technical rows are not statistical replicates |
| `references.json` or YAML | Assembly ID, chromosome sizes/views and SHA256 identities; optional blacklist; exact-grid phasing track when compartments are selected |
| Matrix `.mcool` or `.cool` URI | Raw integer counts, requested existing resolution(s), compatible axes; stored weights and provenance for balanced analyses |
| `comparisons.tsv` | `comparison_id,case_condition,reference_condition,design,inference,modules`; needed for TAD events or condition comparisons, not QC-only/compartment-only runs |

Example comparison row (tab-separated):

```text
comparison_id case_condition reference_condition design inference modules
treated_vs_control treated control ~condition descriptive loops;insulation
```

Use real tab characters or the supplied CSV/TSV writer when creating the file. For TAD-only analysis the condition direction is used; the `modules` field still describes loop/insulation adapters, not a formal TAD model. Mouse and human use the same software, with separate assembly-matched projects and reference assets. Different assemblies/grids cannot be combined in one comparison.

## QC and HiCRep

**Purpose:** QC counts within/between-chromosome contacts and measures distance decay, called P(s). HiCRep measures contact-map similarity between sample pairs; its SCC is a stratum-adjusted correlation score, not a biological P value.

**Select:** `--analysis qc`, `--analysis hicrep`, or `--analysis qc,hicrep`. HiCRep needs its checked isolated environment. HiCRep-only uses raw counts and does not require balancing weights; QC includes balanced P(s) and does.

**Inputs/settings:** QC uses existing `PS_RESOLUTION` (default 10,000 bp), `PIXEL_CHUNKSIZE`, and `INSULATION_IGNORE_DIAGS`. HiCRep uses `HICREP_RESOLUTION=10000`, `HICREP_MAX_DISTANCE_BP=2000000`, and `HICREP_SMOOTHING_BINS=2`. Changing smoothing/distance changes the similarity quantity; inspect sparse/undefined scores.

**Outputs:** `01_qc/contacts/SAMPLE.RESOLUTION.json` contains raw cis/trans totals, exclusions and eligible exposure. `01_qc/ps/SAMPLE.RESOLUTION.tsv` contains distance, raw/balanced sums, averages and eligible-pixel counts. `01_qc/hicrep/by_view.tsv` records view-level SCC evidence; `summary.tsv` contains `sample1,sample2,weighted_scc,valid_views`. The HTML report includes decay curves and, when requested, a similarity heatmap. NA means undefined, not zero similarity.

**Troubleshooting:** missing weights block QC but do not initiate balancing. Too little usable data can produce undefined HiCRep scores. Similarity does not prove that samples are independent biological replicates.

## Compartments

**Purpose/select:** broad chromosome organization and oriented compartment scores; `--analysis compartments`.

**Inputs/settings:** existing `COMPARTMENT_RESOLUTION=100000`, compatible views, stored weights, and an exact-grid registered phasing track. Phasing is a common reference signal used to give compartment scores a consistent sign. `MIN_PHASING_CORRELATION=0.1` rejects weak orientation; `MAX_DENSE_BINS` and the compartment memory reservation bound dense calculations.

**Outputs:** `03_features/compartments/SAMPLE.tsv`, `.eigenvalues.tsv`, `.phasing.json`, and `.bedgraph`. Tables retain components, selected component and oriented scores; QC records orientation correlations. BedGraph is a continuous browser track. Ambiguous orientations remain NA. Compartment comparisons/sign changes are descriptive; no formal compartment FDR is generated.

**Troubleshooting:** a missing/misaligned phasing track or oversized view needs input/view correction. Compartments and saddle require phasing; TAD-only and loop-only analyses do not. Full `practical_2_4kb`/`hsieh_inspired` profiles require compartments and saddle plots in the final report; `--analysis` selections still enable only their prerequisites.

## Insulation, boundaries and TADs

**Purpose/select:** insulation measures contacts across a genomic position; minima can mark separating boundaries. Use `--analysis insulation` for scores/boundaries, or `--analysis tads` for scores, per-replicate/condition domains and descriptive event comparisons.

**Inputs/settings:** weights, existing `INSULATION_RESOLUTION=10000`, compatible views, and for TAD events a comparisons table. Windows default to `50000,100000,250000` bp; `INSULATION_PRIMARY_WINDOW_BP=100000`. `INSULATION_MIN_VALID_PIXEL_FRACTION=0.8` controls coverage eligibility. `INSULATION_BOUNDARY_THRESHOLD="Li"` preserves cooltools' automatic threshold; `Otsu` or a fixed nonnegative number is also supported. Windows must be positive multiples of the grid and include the primary window.

**Outputs:** `03_features/insulation/SAMPLE.tsv` contains scores, valid-pixel fractions, prominence and boundary flags. `SAMPLE.WINDOW.bedgraph` is a browser score track. `03_features/boundaries/SAMPLE.WINDOW.tsv` contains called boundaries. TAD outputs include BED/TSV domains, support/disagreement, event tables, sensitivity results, and report figures. See [TAD parameters, exact files and event meanings](18-tad-analysis.md).

**Troubleshooting:** NA near edges/masks is expected. A missing boundary is not automatically biological disappearance. Higher fixed thresholds generally retain fewer boundaries; larger windows change the scale being measured.

## Loops and quantification

**Purpose/select:** enriched contacts between two genomic positions; `--analysis loops`. `--analysis quantify` also discovers loops if no independent feature registry is supplied; with an external registry it measures that set directly.

**Inputs/settings:** existing `LOOP_RESOLUTION=10000`, weights, compatible views and matching expected contacts calculated from the same masks. Dots discovery automatically enables QC. Defaults retain cis anchors 50,000–2,000,000 bp apart. `LOOP_CONSENSUS_TOLERANCE_BINS=1` matches nearby calls without transitive drift. Optional 5-kb calling requires an existing resolution and explicit reviewed approval; no finer grid is silently built.

**Outputs:** `03_features/loops/SAMPLE.RESOLUTION.bedpe.tsv` holds paired genomic coordinates and available caller statistics; `consensus.tsv` holds common IDs/coordinates and support. `04_quantification/discovered_loops/` and `independent_loops/` contain `counts.tsv`, `exposures.tsv`, `exclusions.tsv`, and `features.json`. Count columns are biological sample IDs. `eligible_cis_exposure` is the full eligible raw-contact denominator, not the sum of selected loop counts. An eligible absent pixel is zero; excluded measurements are NA.

**Troubleshooting/visualization:** a header-only zero-call table can be a valid result. A failed caller is an error. BEDPE-style TSVs contain headers/additional columns; export the required coordinate columns for your browser. Caller enrichment P/q values do not establish differential changes between conditions.

## APA, ADA and saddle individually

O/E means observed contacts divided by expected contacts at the same genomic distance/view. These summaries are descriptive and retained per biological replicate.

| Selection | Inputs/prerequisites | Settings and meaning | Outputs |
|---|---|---|---|
| `--analysis apa` | Loop grids/weights; independent candidates if supplied, otherwise loop discovery/consensus | `LOOP_QUANTIFICATION_RESOLUTIONS`, `APA_FLANK_BP`; average enrichment around anchors | `05_aggregates/apa/RES/SAMPLE.npz` with `mean,sum,valid_feature_count`; per-grid summary, condition means and report PNG/SVG |
| `--analysis ada` | Insulation/consensus boundary-derived domains and `ADA_RESOLUTION` grid/weights | Primary insulation window; domains rescaled to 50 bins; O/E inside the common domain set | `05_aggregates/ada/SAMPLE.npz` with `mean,sum,valid_domain_count`; condition means, report PNG/SVG and `04_quantification/domain_scores/strength.tsv` |
| `--analysis saddle` | Compartments, exact phasing and `COMPARTMENT_RESOLUTION` grid/weights; no loop grid needed | `SADDLE_PARTITION_SOURCE=replicate_rank_mean`; 50 fixed compartment groups; optional registered external/pool source | `05_aggregates/saddle/SAMPLE.npz` with `mean,sum,valid_pixel_count`; condition means, partition TSV/JSON, AA/BB/AB/segregation summary and report PNG/SVG |

`AGGREGATE_MODULES="apa,ada,saddle"` controls aggregates in a configured full pipeline. Full practical/Hsieh-inspired profiles require saddle and its compartment/report prerequisites. Named selection chooses only its requested aggregate consumers. Empty loop/domain feature sets can produce valid arrays containing NA and zero contribution counts, with an explicit report explanation. A saddle without any eligible shared oriented partition fails rather than producing a completed empty plot. Do not treat the number of pileup windows as biological replication.

## Condition comparisons

Inferential loop comparisons automatically include [volcano plots](22-loop-volcano-plots.md) in the enabled report: effect versus biological P value, colours by model FDR, PNG/SVG exports and plotted-point TSVs. Descriptive comparisons explain their absence. Plot settings change display only; they do not alter tests.

**Select:** `--analysis comparisons` uses declared loop/insulation modules; `--analysis loop_comparison` or `--analysis insulation_comparison` restricts the comparison adapters to that module.

**Inputs/settings:** a comparisons table plus relevant measurements. Descriptive loop comparisons use discovered features, or a registered independent set. Formal/preliminary loop comparisons require independently frozen candidates; they do not require same-data discovery. Insulation comparisons use the fixed score grid/primary window. Set `inference` per comparison row, not with `RUN_DIFFERENTIAL` alone. `FDR_THRESHOLD=0.05` and prespecified filtering/calibration settings are documented in [statistics](05-statistics.md).

**Outputs:** `06_comparisons/ID/parameters.json`; model-input TSVs; `loops.tsv` and/or `insulation.tsv`. Descriptive effects are `log2_allocation_ratio` and `case_minus_reference_log2_insulation`; P/FDR are NA. Inferential adapters additionally export uncertainty, model diagnostics, P/FDR and exclusion status. Higher insulation is not automatically a stronger boundary. Formal results estimate relative quantities, not absolute contacts per cell.

**Troubleshooting:** R package availability is only one requirement. Independent units, rank/residual degrees of freedom, independent loop candidates and matching calibration must pass. See [R inspection and acceptance](19-r-environment-and-inference.md). Unsupported formal TAD/compartment analysis cannot be enabled by changing a label.

## Optional pools and depth sensitivity

`--analysis pool` requires `RUN_POOL=true`, an explicit `POOL_MODE`, and valid membership. Prefer existing upstream pools; import verifies members/count sums/weights. Creating new pooled maps requires its separate opt-in for balancing changed counts. Outputs are `02_pools/registry.json`, provenance and explicitly requested maps; pools are discovery inputs, not statistical replicates.

`--analysis depth_sensitivity` requires `RUN_DEPTH_SENSITIVITY=true`. Settings include `DEPTH_SENSITIVITY_MODE`, `DOWNSAMPLING_REPEATS=5`, `RANDOM_SEED`, and retention. Raw-only thinning produces P(s) and optional HiCRep; balanced-feature sensitivity needs explicitly permitted new-map balancing. Results live under `01_qc/depth_sensitivity/`, with seed/depth/status summaries and requested fixed-feature products. Depth matching does not manufacture biological replication. Both optional analyses remain off by default; selecting their names alone does not grant permission to create derived maps.

## Reports, variants and errors

Every selected run retains normalized configuration, input SHA256 identities, software inventories, logs, completion records and measured task timing. `08_reports/index.html` explains selection, parameters, stage status and downloadable products. Keep the complete run folder when downloading it so images and relative links work.

Create a parameter variant beside the original configuration through the same installed command:

```bash
microc2compare -c "$m2c_conf" --variant tads-threshold-05 \
  --set INSULATION_BOUNDARY_THRESHOLD=0.5
```

The command preserves the original file and prints a new sibling configuration with the new run ID and literal parameter override. Use that printed path with `--analysis tads --dry-run`, then `--analysis tads`. Relative paths retain their meaning. No Python/source path is required. The advanced standalone `create_variant.py` helper remains available.

All variants read the same original matrices. Different run IDs preserve earlier results; cross-run scientific product reuse is not automatic. Within one run, parameter/signature checks reuse compatible outputs. Source changes are deliberately conservative and can require broader recomputation. Never edit completion signatures to force reuse.

If a stage fails, inspect `logs/STAGE.SCOPE.log`, correct its cause and resume. Missing weights/resolutions require upstream correction; insufficient memory needs smaller views/workers or a measured reservation. A writable output directory is required; software does not enforce HOME-only paths. [Troubleshooting](10-troubleshooting.md) covers locks, dependencies and other failures. [The master output inventory](09-outputs.md) explains retained products and formats.
