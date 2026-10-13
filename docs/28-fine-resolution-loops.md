# Fine-resolution loops: 2 kb and 4 kb

Release 0.3.2 provides a complete multi-resolution path: Mustache discovery, cross-resolution catalogue, raw counts and exposures, fixed-candidate Chromosight scores, APA and comparison tables at each requested grid. The default analysis chromosome scope is **autosomes**. X/Y are an explicit option. Imported matrices and their weights remain read-only.

This is a practical Hsieh-inspired strategy, not an exact reproduction. The [Hsieh publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC9729117/) supports 2/4-kb loop measurements and the 5-kb–2-Mb EP/PP comparison range. Historical caller thresholds, matching rules and an insulation-calling window have not all been verified. A displayed smoothing window is not evidence of the calling window. Our mm39 data also differ from the publication's mm10 reference.

## 1. Inventory existing matrices

Use your normal Bash login and the accepted shared command. Set the scientific grids and sample/reference pointers once in the project configuration, then:

```bash
microc2compare -c /your/project/config/project.conf --audit-resolutions
```

`availability.tsv` distinguishes missing grids, invalid weights, unverified convergence and metadata-usable inputs; `audit.json` records membership/scope and identities. This metadata audit does not establish biological depth/reproducibility. [Simple project commands](32-simple-project-commands.md) covers setup without versioned script or Python paths. The advanced standalone `audit_resolutions.py` remains available for inspecting exploratory grids outside the configured request.

## 2. Prepare only missing grids when justified

```bash
nohup microc2compare -c /your/project/config/project.conf --prepare-resolutions \
  > /your/project/fine-grid-preparation.log 2>&1 < /dev/null &
```

This explicit operation reuses existing/registered grids and creates only missing exact integer coarsenings in a new project-owned directory. It verifies contact conservation, balances only new matrices with the actual source recipe and checks imported hashes. It prints a new sibling configuration with the augmented sample/provenance and matrix-registry pointers already set. It cannot recover 2-kb data from a 10-kb map; suitable upstream valid pairs are needed then. Failed preparation is not an accepted input set. Biological observations remain separate; no `.hic` regeneration is required.

## 3. Configure a new analysis run

Create a compact configuration with `microc2compare --init /your/project/config/project.conf`, or reuse your existing fine-loop configuration. Set finalized inputs, reference, output and scratch paths. Successful `--prepare-resolutions` already prints a new configuration containing these pointers; no manual assembly is needed. For an advanced manual configuration the equivalent pointers are:

```bash
SAMPLESHEET="/your/project/upstream_regrouping/fine_grids001/samples.tsv"
MATRIX_REGISTRY="/your/project/upstream_regrouping/fine_grids001/matrix_registry.json"
```

The practical profile discovers and measures at 2000 and 4000 bp, uses a 5-kb–2-Mb measurement range, and keeps insulation/TAD/ADA at 10 kb. The example preserves this project's 40-kb minimum TAD length. **The full practical/Hsieh-inspired analysis requires compartments, saddle plots and a final HTML report**, as well as the configured APA/ADA analyses. The example sets `RUN_COMPARTMENTS=true`, `COMPARTMENT_RESOLUTION=100000`, `AGGREGATE_MODULES="apa,ada,saddle"` and `RUN_FINAL_REPORT=true`. Removing saddle or disabling its prerequisites fails configuration validation for a full profile; explicit individual module selections remain available.

Mustache 1.3.3's implementation suppresses the first diagonals. Conservative effective discovery cutoffs are 10 kb at 2 kb and 20 kb at 4 kb, even with a requested 5-kb measurement minimum. The report exposes this limitation. Fixed candidates can be measured closer than the caller can discover; near-diagonal scoring windows can still be unavailable. Do not interpret the lack of a short-range call as absence of promoter–enhancer interaction. Theoretical Micro-C resolution is not guaranteed per-replicate statistical resolution.

### Register the common phasing input before the full run

Compartments need one fixed reference signal to choose a consistent sign: for example, an assembly-matched GC-content track or a reviewed common gene-expression track. Positive oriented scores then follow the chosen signal; it is not a biological replicate or a differential-test input. Reuse a verified existing resource when possible. Do not use each condition's independently signed eigenvector as a common phasing reference.

Register `phasing_track` and `phasing_track_sha256` in a **new copy** of your reference JSON, and point the new configuration's `REFERENCE_REGISTRY` to that copy. This preserves the existing QC run's reference identity. The track is a headerless four-column tab-separated BedGraph (`chrom`, `start`, `end`, numeric value), with zero-based coordinates exactly matching the configured compartment grid, including each chromosome's terminal bin. Do not relabel or interpolate a different genome/grid. Example rows at 100 kb:

```text
chr1	0	100000	0.42
chr1	100000	200000	0.39
```

Those rows illustrate the format; they are not a usable genome-wide reference. Use `microc2compare -c project.conf --register-resources` to calculate/register the actual track's hash. Preflight validates the registered file/hash and coordinate overlap before contact scans. Review per-view `*.phasing.json`; weak orientation remains missing. Saddles require at least 50 commonly usable, consistently oriented bins in a view. If no view qualifies, the saddle stage fails with an actionable explanation instead of reporting an empty plot as a completed analysis.

## 4. Pilot, inspect, then continue

Use the accepted versioned candidate command; see [upgrade instructions](33-simple-project-upgrade.md). The shared alias changes only after server acceptance.

```bash
command="microc2compare"
conf="/your/project/config/project.conf"
"$command" -c "$conf" --analysis qc,hicrep --dry-run
nohup "$command" -c "$conf" --analysis qc,hicrep \
  > "/your/project/fine001-qc.log" 2>&1 < /dev/null &
```

After QC inspect per-resolution contact retention, usable-bin fractions, distance decay and concordance. Then continue with the full configured analysis:

```bash
nohup "$command" -c "$conf" --resume \
  > "/your/project/fine001-full.log" 2>&1 < /dev/null &
```

A changed selection can require additional registrations/QC, which the runner validates. For a faster scientific pilot, register an **explicit** autosomal subset in a separate RUN_ID, such as chr1, using the [scope guide](30-chromosome-scope.md). A parameter variant also needs a new RUN_ID; `--resume` is for continuing its current compatible inputs/settings.

## 5. Separate modules and outputs

| Selection | Inputs beyond common manifests/reference/maps | Main products |
|---|---|---|
| `loops` | Mustache prefix, each discovery grid | `03_features/loops/SAMPLE.RES.bedpe.tsv`, caller JSON, `catalogue.tsv`, `source_membership.tsv` |
| `quantify` | Discovered catalogue or independent registered BEDPE | `04_quantification/{discovered_loops,independent_loops}/RES/{counts,exposures,exclusions,candidate_mapping}.tsv` |
| `chromosight` | Fixed candidate counts, Chromosight prefix | `04_quantification/chromosight/RES/SAMPLE.tsv` and software/kernel JSON |
| `annotation` | GTF and optional cCRE registry | `03_features/annotation/{anchors,loops}.tsv`, `resources.json` |
| `apa` | Fixed mapped candidates and weights | `05_aggregates/apa/RES/*.npz`, summary, condition means, PNG/SVG |
| `tads` / `ada` | 10-kb insulation and configured boundary windows | Existing TAD/event tables; `05_aggregates/ada`, domain strengths and PNG/SVG |
| `compartments,saddle` | Compartment grid and phasing track | Compartment tracks; shared partition, saddle NPZ/summary and PNG/SVG |
| `loop_comparison` | Comparison table; formal candidates/calibration if inference requested | `06_comparisons/ID/RES/loops.tsv`, model inputs/parameters; biological volcano only when appropriate |

For example: `"$command" -c "$conf" --analysis loops,chromosight,apa`. Prerequisites are enabled automatically; unrelated modules are omitted. Use distinct run configurations to preserve parameter variants. Each module produces an HTML report. A full-profile report must contain saddle arrays and PNG/SVG plots for every biological replicate and every equal-replicate condition mean, plus the per-replicate segregation summary. `--analysis saddle` produces these same report products as a separate module; QC-only or loop-only reports do not require unrelated saddle inputs. Download the whole run directory so its relative links and images work offline. The primary-resolution compatibility aliases remain, but `RES` directories are authoritative for multi-resolution results.

Complete-linkage correspondence prevents a chain of nearby calls from becoming one oversized cluster. Original calls and caller significance remain in source membership. A finest observed call represents each cluster. Anchor midpoints map candidates to measurement grids; an exact midpoint on a bin boundary is assigned to the lower bin. Coarse-to-fine projections are labelled and do not create new localisation evidence. Candidates sharing a coarse pixel become **one count/test**, with all aliases retained. Same-bin, masked and otherwise unavailable measurements are NA; an observable zero remains zero.

## 6. Statistical interpretation and optional finer grids

All eligible candidates are measured in every biological observation, whether or not that observation called the loop. Same-data discoveries and condition means support exploratory effects. Biological differential P values require independently frozen candidates, an identifiable design, sufficient independent units and a passing calibration for the resolution/distance/count policy. Caller enrichment q values and Chromosight correlations are not differential FDR. Multiple grids are correlated sensitivity measurements; do not pool their P values.

Keep the four basal observations and the four mESC cell-line units defined by the reviewed sample sheets. Within each mESC unit rep1/rep2 remain technical inputs. Source study and cell type are confounded in this comparison; the descriptive design cannot isolate a causal cell-type effect from study/library differences.

1-kb and sub-kb loop analysis needs `RESOLUTION_ELIGIBILITY_REGISTRY`, binding a current audit and an explicit retention/reproducibility review for every requested fine grid. Sub-kb also requires `ALLOW_SUBKILOBASE_ANALYSIS=true`. The example registry starts unapproved. 2/4-kb metadata acceptance still requires a real-data pilot. No production statistical calibration is supplied.

## Resources and troubleshooting

The example reserves 24 CPUs and 160000 decimal MB on the shared server, with two concurrent matrix tasks and four loop workers each. Mustache processes overlapping bounded dense tiles; its worker estimate must fit `MEMORY_MB_LOOPS` and `MAX_DENSE_BINS`. Dense cost grows roughly with the square of tile width. Lowering resolution or extending maximum distance can therefore increase memory sharply. Threads are configurable but are not guaranteed speedup; annotation/quantification and aggregate branches use less parallelism. Nested BLAS threads stay one. Coordinate other users' budgets; this runner is not a server scheduler.

Missing grid: reuse a true finer grid or rebuild upstream. Failed weights: correct only the relevant new/upstream product; never bypass convergence. Missing Mustache/Chromosight: provision/check the optional loop prefix. Dense memory error: reduce workers/distance, or raise a measured reservation within available memory. Sparse/flat windows: inspect status and retention rather than zero filling. Changed chromosome scope: start a new RUN_ID. Empty loop discoveries are valid negative software outputs; review depth, effective diagonal cutoff and parameters before biological interpretation.
