# Both-anchor gene and cCRE annotation

Annotation is optional and independently selectable: `microc2compare -c project.conf --analysis annotation`. It uses the discovered catalogue or an independently registered candidate set. Technical-library merging remains upstream. Gene names and regulatory labels do not change raw count measurements or define statistical replication.

## Reuse verified resources

Follow the resource conventions in [ATACseq2tracks annotation](https://github.com/MichalGd/ATACseq2tracks/blob/main/docs/16_peak_annotation.md): a registered GTF and an assembly-matched cCRE BED. Reuse readable central files when their genome, release and conversion history are verified; do not download or copy them for every project. On this server `/opt/bioinformatics/ATACseq2tracks_shared/references/mm39/annotation.gtf` and `mm39.ccre.bed` are candidate locations, **not verified server findings**.

Reuse a verified shared annotation registry when available. Otherwise copy `config/annotation.example.json` into your project and supply only the actual paths and source descriptions. Record the GTF release, cCRE release and any liftover source/target assembly, chain and filtering history. Leave missing SHA256 fields empty; do not type checksums by hand. In the reviewed ATAC workflow, mouse cCREs are ENCODE3 mm10 records lifted to mm39; they should not be relabelled native mm39 ENCODE4 data. Human references require their own matching resources. GTF alone is supported: omit cCRE keys when disabled.

Set `ANNOTATION_REGISTRY="annotation.json"`, then run:

```bash
microc2compare -c /your/project/config/project.conf --register-resources
```

Use the printed sibling configuration. This calculates missing hashes and checks existing identities without modifying/copying shared assets. Annotation semantics and conversion history still require review. See [simple project commands](32-simple-project-commands.md).

In `project.conf` set:

```bash
RUN_ANNOTATION="true"
ANNOTATION_REGISTRY="annotation.json"
# Optional overrides; these values are already the internal defaults.
# ANNOTATION_PROMOTER_UPSTREAM="2000"
# ANNOTATION_PROMOTER_DOWNSTREAM="500"
```

The GTF must include gene records with gene IDs and `+`/`-` strand, plus exon records for exon/intron distinctions. GTF one-based inclusive coordinates are converted to zero-based half-open intervals. cCRE BED must have at least chromosome, start, end and ID. Automatic detection recognises the legacy classes `PLS`, `pELS`, `dELS`, `CTCF-only` and `DNase-H3K4me3` in later columns. For another layout or vocabulary, specify `ccre_class_column` as a **zero-based** column index following BED ID. The explicit-column path retains the supplied labels; review their vocabulary rather than assuming it performs biological validation. A row without automatically recognised classes fails if no column is selected. Reference mismatch or checksum drift also fails.

For the reviewed mouse BED, column 10 contains classes (`ccre_class_column: 9`); for the six-column human ENCODE4 file, column 6 contains classes (`ccre_class_column: 5`). Explicit selection retains human `CA-CTCF`, `CA-TF`, `CA`, `CA-H3K4me3` and `TF` labels without recoding them as enhancer-like evidence. [SCREEN class definitions](https://screen.wenglab.org/about) describe these signatures. The installed parser uses the fourth-column BED identifier in `ccre_id`; the reviewed human BED also carries a distinct fifth-column `EH38E...` cCRE accession. Keep that distinction in registry provenance and use the original source file for joins to the additional accession.

## What is recorded

Each distinct anchor receives a stable ID and these annotations:

- Strand-aware promoter overlaps, all overlapping gene IDs/names, exon and gene-body overlaps. The primary context uses promoter, exon, intron, then distal intergenic precedence; all overlap evidence remains available.
- Nearest gene TSS measured from the anchor midpoint. Ties are retained; signed distance follows gene strand. This explicit midpoint convention differs from an interval-distance nearest-feature convention and is labelled in provenance.
- Every overlapping cCRE BED ID/class. A primary match uses greatest overlap, then PLS, pELS, dELS, DNase-H3K4me3, CTCF-only precedence, then stable ID; newer explicit-column classes are retained but unranked in that legacy tie-break order. All alternatives remain available.
- Both anchors on every candidate, with `promoter_promoter`, `promoter_enhancer_like`, `enhancer_like_enhancer_like` or `other_or_unresolved` descriptive categories. pELS/dELS supply enhancer-like evidence. Disabled cCRE annotation stays distinguishable from an enabled file with no overlap.

Products are `03_features/annotation/anchors.tsv`, `loops.tsv`, browser-ready `anchors.bed`, `loops.bedpe`, and `resources.json`. The report explains categories and lists their counts. These annotate original catalogue coordinates, not collapsed coarse measurement pixels; join quantitative tables through `candidate_mapping.tsv`.

cCRE or nearest-gene annotations do not prove an active enhancer, a regulatory target or causality. Cell-type-matched accessibility/H3K27ac evidence can strengthen interpretation, but this release does not silently infer activity from unmatched public resources. Preserve ambiguous promoter overlaps and avoid choosing one gene solely to simplify a figure.
