# Configuration reference

On the shared server, copy project configurations into your own workspace and use the published launcher, which loads the checked software profile internally. Keep results/scratch outside the release; resource presets apply per run. See [shared user setup](16-shared-server-installation.md#3-start-as-an-ordinary-server-user).

The project `.conf` is trusted Bash and is sourced after template defaults. Unknown uppercase assignments in that file are rejected; normalized allowed values are validated and snapshotted. This is a typo check, not a Bash security sandbox. Never source sample rows or an upstream configuration automatically.

Relative configuration paths resolve against the project configuration directory. Manifest paths resolve against their own table directory. CLI `-s` resolves against the invocation directory. Source a resource preset explicitly, then place project overrides afterward. Missing optional paths are empty strings.

```bash
source "/path/to/microc2compare-0.4.0/config/resources_pilot.conf"
PROJECT="experiment"
SAMPLESHEET="samples.csv"
TECHNICAL_MEMBERSHIP="technical_membership.tsv"
COMPARISONS="comparisons.tsv"
REFERENCE_REGISTRY="references.yaml"
RUN_ID="new_variant001" # software settings are supplied by the shared command
```

| Variable | Template default | Meaning |
|---|---|---|
| `PIPELINE_NAME` | `microc2compare` | Workflow identity; keep microc2compare. |
| `SCHEMA_VERSION` | `1` | Configuration schema version; v0.1 accepts 1. |
| `PROJECT` | `project` | Safe project identifier, matching every sample row. |
| `ANALYSIS_SET_ID` | `analysis` | Stable name for the complete common analysis set. |
| `RUN_ID` | `run001` | Stable explicit run identity; change for a new historical run. |
| `REFERENCE_ID` | `reference` | Registered assembly/reference identity; must match every sample. |
| `ANALYSIS_MODE` | `within_assay` | within_assay or explicit descriptive assay_comparison. |
| `SAMPLESHEET` | `samples.csv` | Biological observations CSV/TSV; one consolidated map per unit/condition. |
| `UPSTREAM_INPUTS` | `` | Optional experiment-specific upstream input table for --prepare-inputs; existing finalized manifests can be used directly. |
| `TECHNICAL_MEMBERSHIP` | `technical_membership.tsv` | Technical library/source/deduplication TSV. |
| `COMPARISONS` | `comparisons.tsv` | Explicit contrasts, design, inference mode and module table. |
| `REFERENCE_REGISTRY` | `references.yaml` | Immutable reference JSON/YAML with asset checksums. |
| `INDEPENDENT_FEATURES_REGISTRY` | `` | Frozen external independent loop registry; required for formal/preliminary loops. |
| `OUTDIR` | `../results` | Owned results base; appends project/analysis/run identifiers. |
| `SCRATCH_ROOT` | `../scratch` | Owned temporary base; task-scoped TMPDIR, never imported data. |
| `RANDOM_SEED` | `20261007` | Recorded PCG64 sensitivity/calibration base seed. |
| `PYTHON_BOOTSTRAP` | `python3` | Python 3.11+ executable for dependency-free metadata, and explicit development helpers. |
| `CONDA_EXE` | `conda` | Registered Conda executable; no solve/update during analyses. |
| `ENVIRONMENT_MODE` | `reuse` | reuse immutable checked prefixes; locked rejects drift/missing locks; development explicitly uses local Python/R. |
| `CORE_ENV_PREFIX` | `/opt/conda_envs/microc2compare/envs/core-0.1.0-clone-20261007` | Checked core scientific Conda prefix. |
| `HICREP_ENV_PREFIX` | `/opt/conda_envs/microc2compare/envs/hicrep-0.1.0-shared-20261007` | Checked separate HiCRep Conda prefix, or a compatible checked shared prefix. |
| `DIFF_R_ENV_PREFIX` | `/opt/conda_envs/microc2compare/envs/diff-r` | Checked R/edgeR/limma prefix; needed only for inferential comparisons. |
| `CORE_LOCK` | `` | Explicit tested Linux Conda lock, with matching .inventory.json in locked mode. |
| `HICREP_LOCK` | `` | Explicit HiCRep lock; pip wheel artifacts must also be installed/verified separately. |
| `DIFF_R_LOCK` | `` | Explicit tested R Conda lock, with matching .inventory.json in locked mode. |
| `R_SCRIPT` | `Rscript` | Rscript executable used only in explicit development mode. |
| `RUN_QC` | `true` | Contact summaries and raw/balanced P(s); required by loops/aggregates/sensitivity. |
| `RUN_HICREP` | `true` | Pairwise concordance task at the QC barrier even if balanced QC is disabled; also in enabled depth sensitivity. |
| `RUN_POOL` | `false` | Optional pool discovery maps; default false, no core dependency. |
| `POOL_MODE` | `import` | import existing validated pools, or explicit create for missing requested scopes. |
| `POOLS_MANIFEST` | `` | Optional pool membership/URI/provenance table; mandatory for import mode. |
| `RUN_DEPTH_SENSITIVITY` | `false` | Optional descriptive seeded contact thinning; default false. |
| `DEPTH_SENSITIVITY_MODE` | `raw_only` | raw_only P(s)/HiCRep, or explicit balanced_features fixed-loop/insulation summaries. |
| `BALANCE_NEW_DERIVED_MAPS` | `false` | Allows balancing genuinely changed task-owned optional maps only; alone launches no task. |
| `RUN_COMPARTMENTS` | `true` | Oriented eigenvectors at the configured compartment grid; needs phasing track. |
| `RUN_INSULATION` | `true` | Fixed-grid insulation and per-window boundary calls. |
| `RUN_LOOPS` | `true` | Dots discovery at the primary and approved sensitivity grid. |
| `RUN_AGGREGATES` | `true` | Per-replicate APA/ADA/saddles using enabled feature consumers. |
| `RUN_DIFFERENTIAL` | `true` | Explicit comparison adapters; formal requests retain all statistical gates. |
| `RUN_FINAL_REPORT` | `true` | Standalone HTML, plots, methods and downloadable outputs; required for full practical_2_4kb/hsieh_inspired profiles, including replicate/condition saddle plots. |
| `KEEP_POOLED_MAPS` | `false` | Retain newly created pool maps after consumers; imported pools always stay untouched. |
| `KEEP_DOWNSAMPLED_MAPS` | `false` | Retain explicitly generated sensitivity maps instead of releasing after each seed. |
| `CLEAN_SCRATCH_ON_SUCCESS` | `true` | Remove only the successful task-owned scratch directory. |
| `COMPARTMENT_RESOLUTION` | `100000` | Existing input group for eigenvectors; v0.1 scientific default 100 kb. |
| `INSULATION_RESOLUTION` | `10000` | Existing input fixed score grid; v0.1 default 10 kb. |
| `INSULATION_WINDOWS_BP` | `50000,100000,250000` | Comma-separated positive grid-multiple windows; includes primary. |
| `INSULATION_PRIMARY_WINDOW_BP` | `100000` | Single primary differential insulation family; default 100 kb. |
| `INSULATION_MIN_VALID_PIXEL_FRACTION` | `0.8` | Prespecified eligible/full-diamond pixel coverage threshold. |
| `INSULATION_IGNORE_DIAGS` | `2` | Ignored diagonals for expected/insulation and explicit new-map balancing. |
| `LOOP_RESOLUTION` | `10000` | Primary single-bin loop grid; default 10 kb. |
| `LOOP_SENSITIVITY_RESOLUTIONS` | `` | Empty or 5000 only; requires existing grid and documented approval. |
| `FINE_RESOLUTION_QC_APPROVED` | `false` | Legacy informational flag; does not authorize 1-kb/sub-kb analysis. Use a reviewed eligibility registry. |
| `LOOP_CALLER` | `cooltools_dots` | cooltools_dots on supported >=5-kb discovery grids, or the pinned Mustache adapter for fine grids. |
| `LOOP_MIN_DISTANCE_BP` | `50000` | Inclusive same-view cis lower bound for calls/counts/exposure. |
| `LOOP_MAX_DISTANCE_BP` | `2000000` | Inclusive same-view cis upper bound for calls/counts/exposure. |
| `LOOP_CONSENSUS_TOLERANCE_BINS` | `1` | Complete-linkage anchor tolerance in bins; default one, no transitive chaining. |
| `LOOP_COUNT_WINDOW_BINS` | `1` | Single pixel only in v0.1; values other than one are rejected. |
| `LOOP_MIN_TOTAL_COUNT` | `20` | Prespecified all-sample raw-count sum filter. |
| `LOOP_MIN_MEAN_COUNT_AT_MIN_EXPOSURE` | `5` | Prespecified mean allocation times minimum exposure filter. |
| `LOOP_MIN_FEATURES_FORMAL` | `1000` | Minimum retained independent features before dispersion fitting. |
| `LOOP_NORMALIZATION` | `eligible_cis_exposure` | eligible_cis_exposure only; no TMM/weights/OE in the raw count model. |
| `PS_RESOLUTION` | `10000` | Existing input resolution for P(s) and basic domain registration. |
| `HICREP_RESOLUTION` | `10000` | Existing raw-count resolution for concordance. |
| `HICREP_MAX_DISTANCE_BP` | `2000000` | Maximum SCC distance; views are combined by stratum weights. |
| `HICREP_SMOOTHING_BINS` | `2` | Sparse smoothing half-width; its excluded-bin footprint is masked. |
| `DOWNSAMPLING_REPEATS` | `5` | Inactive unless depth sensitivity is enabled; default five draws. |
| `FDR_THRESHOLD` | `0.05` | Declared model BH threshold; calibration protocol currently fixes nominal 0.05. |
| `FORMAL_MIN_UNITS_PER_CONDITION` | `3` | Minimum independent units for formal requests. |
| `PRELIMINARY_MIN_UNITS_PER_CONDITION` | `2` | Explicit preliminary unit floor, still requires rank/df/calibration. |
| `MIN_RESIDUAL_DF` | `2` | Residual degrees-of-freedom floor; policy requires at least two. |
| `CALIBRATION_REGISTRY` | `../validation/calibration_registry.json` | Hashed passing artifacts keyed to source/environment/module/design/recipe. |
| `ALLOW_UNCALIBRATED_DEVELOPMENT` | `false` | Explicit fixture/method exception; outputs labeled uncalibrated, never release evidence. |
| `WEIGHT_COLUMN` | `weight` | Existing upstream bin column; never written to imported maps. |
| `WEIGHT_CONVENTION` | `multiplicative` | Explicit multiplicative or divisive; divisive conversion is in memory. |
| `LEGACY_WEIGHT_VALIDATION` | `reject` | reject missing convergence, or explicit finite_positive exception with evidence gap. |
| `SADDLE_PARTITION_SOURCE` | `replicate_rank_mean` | replicate_rank_mean, external registered track, or validated all_sample_pool. |
| `SADDLE_REFERENCE_TRACK` | `` | Exact-grid external bedGraph for the explicit external partition mode. |
| `MIN_PHASING_CORRELATION` | `0.1` | Minimum absolute component/track correlation for orientation. A registered common phasing track is required for compartments and saddles; ambiguous views remain missing. |
| `MAX_DENSE_BINS` | `12000` | Bound on dense regional arrays; compartment memory estimate also applies. |
| `PIXEL_CHUNKSIZE` | `1000000` | Raw pixel streaming/merge chunk size; tune with measured memory/I/O. |
| `MIN_FREE_DISK_MB` | `1000` | Minimum free-space floor near outputs/scratch; optional derived capacity needs pilot budgeting. |
| `TOTAL_CPUS` | `24` | Conservative CPU admission budget within actual affinity/cgroup allocation. |
| `MAX_MEMORY_MB` | `160000` | Decimal MB admission budget within actual available physical/cgroup memory. |
| `MAX_PARALLEL_JOBS` | `8` | Total homogeneous batch job cap. |
| `MAX_PARALLEL_SAMPLES` | `2` | Biological-sample/batch cap; conservatively also limits run-level stages. |
| `MAX_PARALLEL_MATRIX` | `2` | Cap on simultaneous matrix-heavy tasks. |
| `MAX_PARALLEL_IO` | `1` | Cap on simultaneous registration/pool/thinning I/O-heavy tasks. |
| `THREADS_QC` | `4` | Conservative qc worker/admission reservation; serial adapters may use less. Nested BLAS threads stay one. |
| `MEMORY_MB_QC` | `16000` | Decimal MB qc memory reservation; admission bound, not enforced RSS. Must fit the global budget. |
| `THREADS_POOL` | `4` | Conservative pool worker/admission reservation; serial adapters may use less. Nested BLAS threads stay one. |
| `MEMORY_MB_POOL` | `32000` | Decimal MB pool memory reservation; admission bound, not enforced RSS. Must fit the global budget. |
| `THREADS_BALANCE` | `8` | Conservative balance worker/admission reservation; serial adapters may use less. Nested BLAS threads stay one. |
| `MEMORY_MB_BALANCE` | `64000` | Decimal MB balance memory reservation; admission bound, not enforced RSS. Must fit the global budget. |
| `THREADS_COMPARTMENTS` | `4` | Conservative compartments worker/admission reservation; serial adapters may use less. Nested BLAS threads stay one. |
| `MEMORY_MB_COMPARTMENTS` | `16000` | Decimal MB compartments memory reservation; admission bound, not enforced RSS. Must fit the global budget. |
| `THREADS_INSULATION` | `4` | Conservative insulation worker/admission reservation; serial adapters may use less. Nested BLAS threads stay one. |
| `MEMORY_MB_INSULATION` | `24000` | Decimal MB insulation memory reservation; admission bound, not enforced RSS. Must fit the global budget. |
| `THREADS_LOOPS` | `8` | Conservative loops worker/admission reservation; serial adapters may use less. Nested BLAS threads stay one. |
| `MEMORY_MB_LOOPS` | `64000` | Decimal MB loops memory reservation; admission bound, not enforced RSS. Must fit the global budget. |
| `THREADS_HICREP` | `2` | Conservative hicrep worker/admission reservation; serial adapters may use less. Nested BLAS threads stay one. |
| `MEMORY_MB_HICREP` | `32000` | Decimal MB hicrep memory reservation; admission bound, not enforced RSS. Must fit the global budget. |
| `THREADS_QUANTIFY` | `2` | Conservative quantify worker/admission reservation; serial adapters may use less. Nested BLAS threads stay one. |
| `MEMORY_MB_QUANTIFY` | `16000` | Decimal MB quantify memory reservation; admission bound, not enforced RSS. Must fit the global budget. |
| `THREADS_AGGREGATES` | `4` | Conservative aggregates worker/admission reservation; serial adapters may use less. Nested BLAS threads stay one. |
| `MEMORY_MB_AGGREGATES` | `32000` | Decimal MB aggregates memory reservation; admission bound, not enforced RSS. Must fit the global budget. |
| `THREADS_DIFFERENTIAL` | `1` | Conservative differential worker/admission reservation; serial adapters may use less. Nested BLAS threads stay one. |
| `MEMORY_MB_DIFFERENTIAL` | `16000` | Decimal MB differential memory reservation; admission bound, not enforced RSS. Must fit the global budget. |
| `THREADS_DEPTH_SENSITIVITY` | `4` | Conservative depth sensitivity worker/admission reservation; serial adapters may use less. Nested BLAS threads stay one. |
| `MEMORY_MB_DEPTH_SENSITIVITY` | `64000` | Decimal MB depth sensitivity memory reservation; admission bound, not enforced RSS. Must fit the global budget. |
| `THREADS_REPORT` | `1` | Conservative report worker/admission reservation; serial adapters may use less. Nested BLAS threads stay one. |
| `MEMORY_MB_REPORT` | `4000` | Decimal MB report memory reservation; admission bound, not enforced RSS. Must fit the global budget. |
| `THREADS_PREFLIGHT` | `1` | Conservative preflight worker/admission reservation; serial adapters may use less. Nested BLAS threads stay one. |
| `MEMORY_MB_PREFLIGHT` | `4000` | Decimal MB preflight memory reservation; admission bound, not enforced RSS. Must fit the global budget. |
| `RUN_TADS` | `true` | Insulation-defined per-replicate/condition domains and descriptive event comparisons; active with insulation. |
| `AGGREGATE_MODULES` | `apa,ada,saddle` | Comma-separated apa,ada,saddle selection; saddle is required by full practical_2_4kb/hsieh_inspired profiles. Explicit --analysis selections run their chosen consumers. |
| `INSULATION_BOUNDARY_THRESHOLD` | `Li` | Li (existing default), Otsu, or a finite nonnegative fixed prominence threshold; numbers use cooltools strength units. |
| `TAD_MIN_LENGTH_BP` | `30000` | Minimum boundary-to-boundary domain length in bp; descriptive policy. |
| `TAD_MAX_LENGTH_BP` | `3000000` | Maximum boundary-to-boundary domain length in bp; descriptive policy. |
| `TAD_BOUNDARY_TOLERANCE_BP` | `10000` | Complete-span consensus and boundary matching tolerance in bp. |
| `TAD_MIN_REPLICATE_FRACTION` | `0.67` | Minimum fraction of biological observations supporting each condition domain/evaluable region; rounded up. |
| `TAD_MIN_RECIPROCAL_OVERLAP` | `0.5` | Minimum overlap divided by each interval length for a one-to-one match. |
| `TAD_SIZE_CHANGE_BP` | `20000` | Minimum absolute length change in bp; growth/shrinkage must also exceed the fractional threshold. |
| `TAD_SIZE_CHANGE_FRACTION` | `0.1` | Minimum absolute fractional length change relative to reference. |
| `TAD_STRENGTH_CHANGE` | `0.25` | Minimum absolute case-minus-reference change in mean endpoint boundary prominence for an otherwise geometrically stable matched domain; descriptive annotation, not intra-domain contact strength. |
| `TAD_SENSITIVITY_THRESHOLDS` | `0,0.2,0.5` | Alternative fixed boundary thresholds for complete TAD call/event sensitivity; same score grid, no matrix copies. |
| `PLOT_LOOP_VOLCANO` | `true` | Report biological differential loop P values versus edgeR log2 allocation effects when inferential results exist; descriptive runs explain absence. |
| `VOLCANO_MIN_ABS_LOG2_EFFECT` | `0` | Nonnegative effect magnitude for highlighting only; does not filter points, change tests, or recompute BH FDR. |
| `VOLCANO_LABEL_TOP_N` | `5` | Label up to this many highlighted loops by model P value; 0 disables labels, maximum 20. All usable tested points remain plotted. |
| `CHROMOSOME_SCOPE` | `autosomes` | autosomes (default), autosomes_plus_sex, or explicit. Effective coordinates are frozen; changing them requires a new RUN_ID. |
| `SEX_CHROMOSOMES` | `` | Exact comma-separated registered sex chromosomes, e.g. chrX or chrX,chrY; empty by default. |
| `ANALYSIS_CHROMOSOMES` | `` | Exact chromosome names for explicit subset scope; empty for other policies. |
| `ANALYSIS_PROFILE` | `existing` | existing, practical_2_4kb, hsieh_inspired, hsieh_8grid or short_range_4kb. Fine-loop full profiles require compartments/saddle/report; module selections remain independent. Adapted presets do not claim exact historical reproduction. |
| `MATRIX_REGISTRY` | `` | Optional resolution-specific map path overrides, preserving one biological observation per manifest row. |
| `RESOLUTION_ELIGIBILITY_REGISTRY` | `` | Reviewed current metadata/retention/reproducibility evidence required for 1-kb/sub-kb loop analysis. |
| `HSIEH_METHODS_REGISTRY` | `` | Primary-source verified parameters, source locations and explicit differences from the publication. |
| `ALLOW_SUBKILOBASE_ANALYSIS` | `false` | Explicit optional sub-kb permission; does not replace eligibility review or statistical calibration. |
| `LOOP_DISCOVERY_RESOLUTIONS` | `` | Comma-separated discovery grids; empty uses primary plus legacy sensitivity grids. |
| `LOOP_QUANTIFICATION_RESOLUTIONS` | `` | Comma-separated complete count/scoring/APA/comparison grids; empty uses primary. |
| `LOOP_ENV_PREFIX` | `` | Optional separately checked Mustache/Chromosight prefix; empty checks the core for required tools without updating it. |
| `LOOP_LOCK` | `` | Optional exact caller-prefix Conda lock plus reviewed pip requirements/inventory. |
| `MUSTACHE_P_THRESHOLD` | `0.05` | Mustache enrichment significance threshold. Caller significance is not biological differential FDR. |
| `MUSTACHE_SPARSITY_THRESHOLD` | `0.88` | Pinned Mustache local sparsity threshold; practical default 0.88, not a verified historical Hsieh setting. |
| `MUSTACHE_SIGMA_ZERO` | `1.6` | Starting Gaussian scale in bins; pinned practical default 1.6. |
| `MUSTACHE_OCTAVES` | `2` | Number of increasing Gaussian scale octaves; positive integer. |
| `MUSTACHE_ITERATIONS` | `10` | Pinned 1.3.3 library hardcodes ten scale iterations; other values fail. |
| `CROSS_RESOLUTION_TOLERANCE_BINS` | `1.0` | Complete-linkage anchor matching tolerance measured on the coarser grid; original calls and aliases retained. |
| `RUN_CHROMOSIGHT` | `false` | Measure the fixed-candidate kernel correlation in every biological replicate; descriptive scores only. |
| `CHROMOSIGHT_MAX_MISSING_FRACTION` | `0.1` | Maximum missing fraction within the kernel window; failing windows remain NA with a reason. |
| `RUN_ANNOTATION` | `false` | Optional both-anchor gene/cCRE annotation, preserving ambiguous overlaps. |
| `ANNOTATION_REGISTRY` | `` | GTF/cCRE paths, SHA256 identities, reference/source releases and any assembly-conversion history. |
| `ANNOTATION_PROMOTER_UPSTREAM` | `2000` | Strand-aware promoter extension upstream of gene TSS in bp. |
| `ANNOTATION_PROMOTER_DOWNSTREAM` | `500` | Strand-aware promoter extension downstream of gene TSS in bp. |
| `ADA_RESOLUTION` | `10000` | ADA grid independent of fine loop grids; default 10000 bp. |
| `APA_FLANK_BP` | `100000` | Genomic flank around each loop anchor; divisible by every measured APA grid. |
| `THREADS_TADS` | `2` | Conservative tads worker/admission reservation; serial adapters may use less. Nested BLAS threads stay one. |
| `MEMORY_MB_TADS` | `8000` | Decimal MB tads memory reservation; admission bound, not enforced RSS. Must fit the global budget. |
| `RUN_CHROMOSIGHT_DISCOVERY` | `false` | Independent CPU loop-template discovery; separate from fixed-candidate Chromosight scoring. |
| `CHROMOSIGHT_DISCOVERY_RESOLUTIONS` | `` | Comma-separated native discovery grids; empty inherits loop discovery grids. |
| `CHROMOSIGHT_DISCOVERY_MIN_DISTANCE_BP` | `0` | Discovery lower bound; zero inherits LOOP_MIN_DISTANCE_BP. Must remain within the common measurement universe and upstream retained contacts. |
| `CHROMOSIGHT_DISCOVERY_MAX_DISTANCE_BP` | `0` | Discovery upper bound; zero inherits LOOP_MAX_DISTANCE_BP. Cannot exceed the common measurement universe. |
| `CHROMOSIGHT_DISCOVERY_MIN_PEARSON` | `0.3` | Minimum template correlation for independent candidates; descriptive caller scores, not biological P values. |
| `THREADS_CHROMOSIGHT` | `4` | Conservative chromosight worker/admission reservation; serial adapters may use less. Nested BLAS threads stay one. |
| `MEMORY_MB_CHROMOSIGHT` | `16000` | Decimal MB chromosight memory reservation; admission bound, not enforced RSS. Must fit the global budget. |
| `RUN_STRIPES` | `false` | Optional CPU StripePy discovery and condition support on common eligible chromosome segments. |
| `STRIPE_RESOLUTIONS` | `4000` | Comma-separated authentic stripe grids; default 4000 bp. |
| `STRIPE_ENV_PREFIX` | `` | Explicit separate StripePy prefix; never update the core prefix to add the caller. |
| `STRIPE_LOCK` | `` | Optional exact StripePy Conda lock with reviewed hashed pip requirements. |
| `THREADS_STRIPES` | `4` | Conservative stripes worker/admission reservation; serial adapters may use less. Nested BLAS threads stay one. |
| `MEMORY_MB_STRIPES` | `24000` | Decimal MB stripes memory reservation; admission bound, not enforced RSS. Must fit the global budget. |
| `STRIPE_GENOMIC_BELT_BP` | `5000000` | Maximum native caller stripe/contact belt in bp. |
| `STRIPE_MAX_WIDTH_BP` | `100000` | Maximum native stripe width in bp. |
| `STRIPE_MIN_RELATIVE_CHANGE` | `0` | Native percentage descriptor threshold; retained stripes also require inner mean above outer mean. |
| `STRIPE_NORMALIZATION` | `none` | Only none is supported by this tested raw-contact adapter. |
| `STRIPE_MASK_POLICY` | `common_segments` | common_segments: split views at excluded bins; never replace missing biological measurements with zeros. |
| `RUN_BROWSER_EXPORT` | `false` | Export a dated unified browser bundle after enabled scientific stages; do not rerun measurements. |
| `EXPORT_UCSC` | `true` | Write custom-track headers and native BED/interact/bedGraph companions. |
| `EXPORT_HIGLASS` | `false` | Run available bigWig/beddb/bed2ddb converters; report missing converters explicitly. Does not register or render on a server. |
| `BROWSER_OUTPUT_DIR` | `` | Optional new explicit bundle directory; empty uses a dated child of run/browser_tracks. |
| `HIGLASS_CHROMSIZES` | `` | Optional default coordinate profile, validated against frozen reference evidence; native matrix profiles remain separate. |
| `THREADS_BROWSER_EXPORT` | `1` | Conservative browser export worker/admission reservation; serial adapters may use less. Nested BLAS threads stay one. |
| `MEMORY_MB_BROWSER_EXPORT` | `4000` | Decimal MB browser export memory reservation; admission bound, not enforced RSS. Must fit the global budget. |
| `AGGREGATE_WORKER_MEMORY_MB` | `0` | Decimal MB reserved per observation worker; zero chooses a bounded automatic estimate. |
| `AGGREGATE_CHECKPOINTS` | `true` | Reuse only individually committed, identity- and checksum-validated observation/module/grid pileups. |
| `EXPECTED_CACHE` | `true` | Reuse run-owned expected curves only with matching input, weights, masks, views, parameters, implementation and environment identities. |
| `LIMIT_MATRIX_READERS` | `false` | Optional conservative inner/outer reader cap using MAX_PARALLEL_IO; disabled by default. |
| `QUANTIFY_WORKER_MEMORY_MB` | `0` | Decimal MB reserved per raw-count worker; zero chooses a bounded automatic estimate. |
| `BROWSER_ENV_PREFIX` | `` | Optional separate converter environment; preserve existing scientific prefixes. |
| `BROWSER_LOCK` | `` | Optional exact converter Conda lock with reviewed hashed pip requirements. |
| `QC_WORKER_MEMORY_MB` | `0` | Decimal MB reserved per QC grid worker; zero estimates native bins and real pixel chunk without clipping to fit. |

Optional flags do not imply new-map permission. Invalid required weights/resolutions fail with upstream remediation. Statistical policy floors cannot be weakened below three formal/two preliminary units and two residual df. Tuning scientific settings changes completion/calibration keys and requires a reviewed pilot; no setting proves adequate 5-kb coverage automatically.
