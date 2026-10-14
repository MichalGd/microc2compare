# Compare common chromosomes, including X and Y

For release 0.3.0, [explicit chromosome selection](30-chromosome-scope.md) supersedes implicit inclusion: autosomes are the default and X/Y require opt-in. Historical 0.2.x runs below retain their original versioned semantics and should be resumed with their original launcher.

Version 0.2.6 supports imported maps with different native chromosome order or additional scaffolds. It compares exact registered genomic coordinates and intersects usable bins across biological observations. Imported counts and weights remain unchanged; normal analysis does not export new contact matrices.

## Define the analysis regions

The implementation is genome-independent. Mouse can select chr1–chr19 plus X/Y (21 chromosomes); **human can select chr1–chr22 plus X/Y (24 chromosomes)**. Register the actual assembly's chromosome sizes and matching views. `hg38`/`GRCh38` are supported reference aliases in the input helper; do not mix hg19 and hg38. Other references are usable with their own verified registry and exact names/coordinates. There is no hard-coded autosome-count limit in the workflow.

The existing reference registry's `views` BED4 selects chromosomes and regions. Include `chrX` and `chrY` as explicit rows to retain them. Mouse example (the lengths below are placeholders; obtain actual lengths from the verified mm39 reference):

```text
chr1    0    <verified_chr1_length>    chr1
...
chr19   0    <verified_chr19_length>   chr19
chrX    0    <verified_chrX_length>    chrX
chrY    0    <verified_chrY_length>    chrY
```

Use tabs, no header, and real integer lengths. Register the sizes/views paths and their SHA256 checksums in `references.json`; set `REFERENCE_REGISTRY` in the project configuration. A ready command for the basal/mESC project is in its installation guide. There is no extra chromosome-selection command-line flag and no automatic dropping of a requested chromosome.

For a human project, generate the analogous 24 BED4 rows from your verified reference file:

```bash
awk 'BEGIN {OFS="\t"} $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {print $1,0,$2,$1}' \
  /path/to/verified/hg38.chrom.sizes > hg38.autosomes_XY.views.bed
wc -l hg38.autosomes_XY.views.bed
# Expect 24 rows; register this file and its SHA256 in the human reference registry.
```

Use an explicit reviewed path, write a new view file, and use the exact chromosome naming found in the human matrices. These human commands are independent of the mouse-specific basal/mESC project.

Each selected chromosome must exist in **every** biological map with the registered length, resolution and exact bin boundaries. Missing selected chromosomes, different assemblies/coordinates, noninteger counts, failed weights, changed MAPQ/distance filters or incompatible numerical balancing recipes still fail preflight. A chromosome with no common eligible bins remains represented as unavailable; it is not padded with zero observations. An entirely empty common domain fails.

Chromosome-size registry entries that have no view are not requested for analysis. Native sequences outside the views may remain in the files. The workflow maps the registered coordinates onto each file's native IDs before fetching matrices, streaming counts or annotating loop pixels. Insulation and compartment tables are exported on the common registered order, so downstream TAD and condition comparisons use the same grid.

## Weights and scientific interpretation

Masking extra scaffolds does **not** recalculate weights or make upstream balancing scopes identical. Keep actual `filters.chromosome_policy` in each provenance file. This is the only filter field allowed to differ; other contact filters and numerical balancing recipes must still agree. Review this difference during the descriptive real-data pilot. Do not label global normalization as matched solely because numerical parameters match.

X/Y analyses need the sample's sex or karyotype/copy-number information and sufficient coverage. A female or poorly covered sample can leave few or no common usable Y bins. Sparse, excluded or undefined measurements remain NA. Record coverage and avoid interpreting insufficient measurements as disappearing loops/TADs. Inclusion of X/Y does not change the unpaired biological-replicate design.

## Outputs and operation

```bash
microc2compare -c project.conf --dry-run
microc2compare -c project.conf --analysis qc,hicrep
microc2compare -c project.conf --resume
```

`00_metadata/masks/RESOLUTION.scope.json` records selected chromosomes, common usable-bin counts per chromosome, native chromosome lists, excluded sequences, actual upstream policies and balancing recipes. The HTML report and `08_reports/methods.md` explain the scope. `run_summary.json` includes these records.

QC `raw_total`, `raw_cis` and `raw_trans` describe the **entire native input matrix**, so they may include extra scaffolds. Eligible cis exposure and P(s) distance-decay tables use the **common registered views**. `excluded_bin_fraction` is measured within the registered views; `excluded_native_bin_fraction` additionally reflects excluded native sequences. These scope distinctions prevent extra scaffolds from silently changing the differential exposure denominator.

The optional `POOL_MODE=create` pools raw counts on the registered coordinate grid, with chunked sorting and summation. It requires explicit permission through the existing `RUN_POOL`/`BALANCE_NEW_DERIVED_MAPS` settings; these new task-owned pools need new balancing. Imported pools still require exact full native axes shared with all members for complete raw-sum verification. If members have different native scopes, request a new coordinate-aligned pool instead. Pool creation remains disabled for the initial basal/mESC analysis.

Optional depth sensitivity writes only explicitly requested task-owned maps on the common grid. No imported matrix is repaired, balanced or copied by either path. Common-scope support does not relax calibration or independent-feature requirements for formal differential inference.

## Why identical settings can produce different scopes

Configuration variables only take effect if the executing script implements them. The older installed microc2tracks merge helper used for the mESC regrouping did not implement the newer canonical-chromosome filter; it built bins from the full supplied chromosome-size file. The basal maps contain 22 canonical chromosomes, while the new mESC maps contain those 22 plus 39 scaffolds. Their shared chromosome lengths match. This is consistent with different helper behavior, despite the same declared reference/settings, and is not evidence that the shared chromosomes were aligned to a different genome.

An `unknown` coarse-group assembly label is missing metadata, not a verified assembly identity. The input helper now requests provenance/coordinate review for it instead of treating the word `unknown` as a different genome. A known conflicting assembly label still blocks discovery.

## Validation and upgrading

The release includes numerical tests for reordered axes, extra scaffolds, X/Y, common masks, raw counts/exposures, expected tables, raw thinning and optional pool count conservation. The Linux acceptance script runs actual cooltools/HiCRep and TAD stages on new synthetic fixtures with differing scopes. These are software checks, not statistical calibration or validation of the real biological comparison.

Install as a new immutable source release, reuse existing software environments/profile, then run:

```bash
bash "$m2c_src/scripts/validate_common_chromosome_demo.sh" \
  /opt/conda_envs/microc2compare/site/shared-r-ql-fix-20261008.conf \
  "$HOME/microc2compare-work/validation/common-chromosomes-026-$(date +%Y%m%d-%H%M%S)"
```

Use a new run ID for the real-data pilot. Source signatures invalidate incompatible previous results; do not manually copy completion markers. The project-specific guide provides transfer, installation, disconnection-safe acceptance and activation commands.
