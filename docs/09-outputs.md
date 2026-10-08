# Outputs and interpretation

Every run is scoped by `OUTDIR/PROJECT/ANALYSIS_SET_ID/RUN_ID`. Imported contact maps are registered by path/group and SHA256 rather than copied into this tree.

Use [the module catalogue](17-analysis-modules.md) for launch commands, required inputs and parameter changes. [The TAD guide](18-tad-analysis.md#exact-outputs) lists domain/event files and column meanings. Only enabled/selected modules emit scientific products.

| Path | Product and interpretation |
|---|---|
| `00_metadata/config.json` | Normalized configuration, with resolved paths |
| `00_metadata/samples.json` | Biological samples, comparisons and technical membership |
| `00_metadata/inputs.json` | Registered file SHA256, size and mtime identities |
| `00_metadata/environments.json` | Actual package/runtime inventories |
| `00_metadata/preflight.json` | Validated software/input status and allocation checks |
| `00_metadata/masks/*.npz` | Common mask/view labels and per-sample inherited weights, no contact matrix |
| `00_metadata/masks/*.bins.tsv`, `*.views.tsv` | Exact bin axes and analysis views |
| `00_metadata/status/*.json` | Disabled stage states |
| `01_qc/contacts/*.json` | Raw cis/trans totals, eligible loop exposure and excluded-bin fraction |
| `01_qc/ps/*.tsv` | View/distance raw/balanced sums, eligible denominator and averages |
| `01_qc/hicrep/by_view.tsv`, `summary.tsv` | Pair/view SCC, weights, sparse status and combined pair concordance |
| `02_pools/registry.json`, pool provenance | Optional member/URI registration and explicit ownership |
| `02_pools/*/*.cool` | Explicitly created pooled maps, removed after consumers unless retained |
| `03_features/compartments/*` | Oriented/eigenvector TSV, eigenvalues, phasing QC and bedGraph |
| `03_features/insulation/*` | Fixed-grid scores, coverage counts/fractions, prominence, boundary flags and per-window bedGraph |
| `03_features/tads/samples/*`, `conditions/*` | Per-replicate and condition-supported TAD TSV/BED files; unsupported tables preserve disagreement |
| `03_features/tads/parameters.json`, `sensitivity*` | TAD settings/condition-ID mapping and complete fixed-threshold sensitivity evidence |
| `06_comparisons/ID/tad_events.WINDOW.*` | Descriptive event groups, browser intervals and counts: merge/split/grow/shrink/gain/loss/shift/strength/uncertainty |
| `03_features/boundaries/*` | Per-sample/window boundary calls |
| `03_features/loops/*.bedpe.tsv` | Dots calls at primary/explicit sensitivity resolution; discovery significance only |
| `03_features/loops/consensus.tsv`, `boundaries/consensus.*.tsv`, `boundary_derived_domains/consensus.tsv` | Deterministic consensus loops/boundaries and adjacent domains |
| `04_quantification/*/counts.tsv` | Raw biological-sample counts on independent/discovered loop universes |
| `04_quantification/*/exposures.tsv` | One full eligible-cis raw denominator per biological sample |
| `04_quantification/*/exclusions.tsv` | Feature/sample exclusion reasons; masked values remain NA |
| `05_aggregates/apa`, `ada`, `saddle` | Per-replicate O/E means/count arrays and descriptive summaries |
| `05_aggregates/saddle/partition.tsv`, `partition.json` | Fixed view/group membership, hash and too-few-bin view flags |
| `06_comparisons/ID/*` | Model inputs/parameters and module-specific effect/P/FDR/status tables |
| `01_qc/depth_sensitivity/*` | Optional seeded depth/P(s)/HiCRep/fixed-feature summaries and requested retained maps |
| `08_reports/index.html`, figures, methods | Navigable report with downloads and scientific limitations |
| `08_reports/figures/volcano.ID.png`, `.svg` | Differential loop effect versus biological P value; colours use model FDR |
| `08_reports/volcano/ID.points.tsv`, `status.json` | Original plotted statistics/display annotations and availability/filtering counts; no plot for descriptive comparisons |
| `completion/SCOPE/*.done.json` | Validated success/success-empty signatures and product hashes |
| `logs/*`, `benchmarks/*`, `run_summary.json` | Per-task diagnostics, measured wall/RSS evidence and run state |

Exact emitted filenames are defined by the corresponding stage adapter; a disabled module has no scientific placeholder files. Schemas in `schemas/` document manifest/provenance contracts. TSV headers describe measured quantities explicitly. `.npz` aggregate files contain `mean` and contribution counts; missing aggregate pixels remain NaN. BedGraph omits unavailable values and never zero-fills excluded bins; an entirely unavailable track has a comment indicating no finite measurements.

## Worked interpretation

Inferential `06_comparisons/ID/loops.tsv` tables include `dispersion` (the NB dispersion stored by edgeR), `working_NB_dispersion` (the NB dispersion used in its actual fit), `QL_dispersion` (moderated quasi-likelihood dispersion), and `QL_dispersion_field` (the source field name). A shared scalar NB value is correctly repeated across retained features. `SE_log2_approx_QL` and the approximate confidence intervals use the working dispersion and remain explicitly approximate; `Pearson_residual_*` are NB residuals using that working dispersion, not studentized or QL-scaled residuals. These columns are diagnostic, not additional significance tests. Descriptive comparisons have no fitted-model dispersion or uncertainty.

If a loop's raw count rises from 100 to 180 but its sample exposure doubles, its allocation falls from `100/L` to `180/(2L)`: the log2 allocation ratio is approximately −0.15. Raw counts rose, but the configured estimand fell. A dot caller's enrichment P value says nothing about biological differential FDR. Same-data consensus counts report descriptive allocation effects with biological P/FDR missing.

A positive `case_minus_reference_log2_insulation` means a higher relative insulation score in the case condition at that fixed bin; it does not automatically mean a stronger boundary. Prominence exists only at called minima. A significant bin can be annotated as near a boundary afterward; the testing family still includes all eligible primary-grid bins.

A compartment sign switch is descriptive and meaningful only where both components were confidently oriented using the same track. Saddle AA/BB/AB ratios and APA/ADA strengths summarize matched O/E with recorded denominators; compare per-replicate behavior rather than treating pool/seed products as biological observations.

Review `inference_status`, `validation_status`, eligibility reasons and report stage states before using an effect. `uncalibrated` labels can never be converted into production evidence merely because the command exited successfully. Pointwise confidence intervals have no simultaneous-coverage guarantee; loop intervals additionally carry an approximate-method label.
