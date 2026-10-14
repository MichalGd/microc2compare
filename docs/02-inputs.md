# Inputs and biological identity

Every manifest is parsed as CSV or TSV with proper quoting, UTF-8 BOM handling and CRLF support. Headers are required; unknown optional columns are preserved. A row is never evaluated as shell code. Paths in a manifest resolve against that manifest's directory. Paths assigned in the project `.conf` resolve against the configuration directory; CLI `-s` resolves against the invocation directory. Use Linux paths on the target server; Windows paths are local development inputs only.

## samples.csv

| Field | Meaning |
|---|---|
| `sample_id` | Unique stable identifier using letters, digits, `_`, `-`, `.`; reserved names `mask`, `labels`, `pool__*` are rejected |
| `project` | Must match `PROJECT` |
| `condition` | Experimental condition/time-point label |
| `assay` | `microc` or `hic`; mixed assays require `ANALYSIS_MODE=assay_comparison` and descriptive comparisons |
| `biological_replicate` | Biological replicate label; a library/lane is not a replicate |
| `subject_id` | Stable donor/animal identity when known; for cultures with no shared donor identity, use the independently established biological unit identity |
| `experimental_unit_id`, `culture_id` | Optional aliases for an explicitly established independently grown culture within a subject; if both are filled they must agree |
| `reference_id` | Must match the configured registered reference |
| `matrix_path` | Existing `.norm.mcool`, `.cool`, or resolution-specific Cooler URI; one consolidated biological map per row |
| `provenance_path` | JSON sidecar documenting upstream processing and balancing |
| `batch` | Optional batch column, required when included in the design |
| additional covariates | Optional columns used by supported additive design formulas |
| `pairs_path`, `pair_stats_path` | Optional provenance/contact files; no contact-level operation is implemented in v0.1 |

Multiple cultures from one donor in one condition can be registered for descriptive analysis when every row has a distinct explicit `experimental_unit_id` or `culture_id`. Keep their actual shared donor in `subject_id`; a new culture ID does not create a new donor. Each culture ID has one subject owner and occurs at most once per condition. A culture split across conditions may retain its ID. Formal/preliminary comparisons with multiple observations for one subject within a condition fail because the current R adapters do not fit a nested donor/culture model. Repeated subjects across conditions require explicit subject terms for inference. Multiple independent libraries from the same culture belong in one already consolidated sample. Duplicate matrix URIs, source IDs and libraries assigned to multiple samples fail validation. Distinct matrices or culture labels cannot prove biological independence: the experimental and membership/provenance evidence must be reviewed.

## technical_membership.tsv

Required fields are `sample_id`, `library_id`, `source_id`, `source_kind`, `dedup_group_id`, `dedup_scope`, `source_checksum`. `source_kind` is `lane`, `resequencing`, or `independent_library`. Each constituent source occurs once. Checksums identify the upstream source artifact; they are evidence supplied by the upstream manifest, not substituted with an unrelated matrix hash in real projects.

Lane/resequencing rows for the same library require `dedup_scope=joint_library` and the same dedup group. Separately deduplicating lanes and then merging matrices can preserve cross-lane duplicates. Correct that upstream; matrix merging cannot repair it. Independent libraries from one biological unit may be separately deduplicated, but their biological sample identity stays consolidated. Unknown dedup scope blocks formal inference.

## comparisons.tsv

Required columns: `comparison_id`, `case_condition`, `reference_condition`, `design`, `inference`, `modules`. Effects are case minus reference, or case/reference on a log2 allocation scale. `inference` must explicitly be `descriptive`, `preliminary`, or `formal`; invalid formal requests fail without switching modes. `modules` is `loops`, `insulation`, or the semicolon-separated combination. Header-only comparisons are allowed when differential analysis is disabled.

Supported formulas are additive `~condition`, `~batch+condition`, or `~subject_id+condition`. Subject, batch, assay and condition are factors even when their labels contain only digits. Other numeric covariates are continuous; other nonnumeric covariates are factors. Interaction, executable expressions and arbitrary R formulas are rejected. Full rank and at least two residual degrees of freedom are required for inference. See [statistics](05-statistics.md).

## Reference registration

Use one JSON/YAML object with `reference_id`, `chromsizes`, `chromsizes_sha256`, `views`, `views_sha256`. Optional `blacklist`, `phasing_track` assets each require the corresponding SHA256. Chromosome sizes are two tab-separated columns without headers; views are BED4 (`chrom start end name`) without headers. View names must be unique, intervals nonoverlapping, within chromosomes, and complete-bin aligned at requested grids. Every view chromosome must exist in every biological matrix with the registered length and exact bin boundaries. Native chromosome order and extra unselected sequences may differ: coordinates are mapped to each file's native bin IDs. The reference order defines the exported analysis grid. See [common chromosomes and X/Y](27-common-chromosomes.md).

The phasing asset is an immutable, assembly-matched bedGraph (`chrom start end value`) on the exact compartment grid. Compartments require it. Blacklist BED3 excludes any overlapping bin; a common bin mask also intersects all relevant samples' finite positive inherited weights. No downstream mask changes the upstream balancing solution. Chromosome arms can be registered as separate views to avoid crossing centromeres. Bins outside all views are excluded.

## Upstream provenance.json

Required keys: `reference_id`, `count_units="unique_contacts"`, `filters`, `balancing`, `dedup_scope`, `library_ids`, `source_ids`. `filters` declares `mapq`, `min_distance_bp`, `chromosome_policy`; record other consequential filtering keys as well. Compatible samples must have identical contact filters and per-resolution numerical balancing recipes. The chromosome policy may differ when exact common view coordinates pass validation; retain each actual policy in provenance and review its effect on upstream balancing. Never replace it with a fabricated common policy. The biological membership must match the sidecar library/source IDs. Record the upstream commit, commands, software, assembly source and actual input identities when available.

`balancing` is keyed by resolution as strings. Each entry declares `weight_column`, `convention` (`multiplicative`/`divisive`), `converged` and `recipe`. Convergence must be true. Missing legacy convergence can use the explicit `LEGACY_WEIGHT_VALIDATION=finite_positive` exception, with a visible evidence gap; explicitly failed convergence never passes. Infinite/nonpositive weights fail, NaNs represent excluded bins. Divisive weights are converted in memory only. A `.norm.mcool` stores raw integer pixels and weights together; import registers the existing groups, not separate raw/balanced exports.

## Independent feature registry

Formal loops require a JSON/YAML registry for a frozen BEDPE TSV with `feature_set_id`, `path`, `sha256`, `reference_id`, `resolution`, `coordinate_convention="0-based-half-open"`, `source`, `discovery_samples`, `independent_of_analysis_samples=true`, and inspectable `independence_evidence`. BEDPE needs `chrom1 start1 end1 chrom2 start2 end2` headers. Anchors must be one declared-resolution bin, exactly aligned. Duplicate pixels are canonicalized deterministically; a supplied feature ID is replaced by the stable coordinate-derived ID.

External means independent discovery with evidence, not merely a file outside this checkout. Discovery sample IDs cannot overlap analysis sample IDs. Same-data consensus or a condition-blind pool remains descriptive.

## Optional pools.tsv

Required columns: `pool_id`, `pool_scope` (`condition`/`all_samples`), `condition` (empty for all samples), `member_sample_ids` separated by semicolons, `reference_id`, `matrix_path`, `provenance_path`. Scope membership must exactly match the registered biological set; pools never enter biological models. Pool provenance additionally requires `summation="raw_integer_sum"`, exact `member_sample_ids`, and `member_hashes` keyed by sample ID with matrix SHA256 values. Axes, filters, weight metadata, raw totals and **every raw pixel sum** are checked against members. This streamed validation can be expensive but writes no replacement matrix.

`POOL_MODE=import` reuses the registered maps. Explicit `create` may fill missing condition/all-sample scopes at requested resolutions; a one-member scope aliases the sample. New balanced discovery requires its separate balancing switch. Imported pools are never repaired or removed by cleanup.

Use [the upstream preparation helper](13-upstream-handoff.md) to relate your actual microc2tracks sheet/run to these tables. It suggests consolidated maps and exports draft inputs without inventing missing biological, deduplication or filtering evidence. Review its inventory and finish the drafts before normal preflight.

See [editable examples](../config/samples.example.csv) and the self-contained [fixture generator](../scripts/generate_demo.py). Example paths/hashes must be replaced with real upstream evidence; they are not valid project registrations.
