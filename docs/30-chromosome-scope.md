# Autosome default and optional sex chromosomes

Release 0.3.0 defaults to `CHROMOSOME_SCOPE="autosomes"`. Built-in UCSC names select mouse mm10/mm39 chr1–19 and human hg19/hg38/GRCh37/GRCh38 chr1–22. X/Y, mitochondrial sequences and other scaffolds are excluded from analysis unless explicitly selected. Original files, chromosome order and stored weights are preserved.

For other naming styles (`1` instead of `chr1`) or custom references, register exact chromosome classes in `references.json`:

```json
"chromosome_classes": {"autosomes": ["1", "2"], "sex": ["X", "Y"]}
```

This is a syntax example, not a complete mouse/human chromosome list. Classes must be disjoint, unique and present in registered chromosome sizes. Selected chromosomes also need valid registered views and matching lengths/bin coordinates in **every** biological map and requested resolution. Missing requested chromosomes fail; the workflow never silently shortens the requested set.

To include X only:

```bash
CHROMOSOME_SCOPE="autosomes_plus_sex"
SEX_CHROMOSOMES="chrX"
```

To include both, set `SEX_CHROMOSOMES="chrX,chrY"`. Consider sample sex, copy number and sparse coverage. Missing measurements stay NA, rather than zero or biological disappearance. For a small autosomal pilot:

```bash
CHROMOSOME_SCOPE="explicit"
ANALYSIS_CHROMOSOMES="chr1,chr2"
SEX_CHROMOSOMES=""
RUN_ID="fine_chr1_chr2_pilot001"
```

Registered sizes/views can retain a broader inventory: selection is applied without copying or rebalancing maps. Changing scope, chromosome lengths or selected view intervals requires a **new RUN_ID**. The effective scope is frozen before saved configuration is changed; changing scope on resume is refused. Older releases did not freeze this new policy, so resume their historical runs with the original versioned command. Use a new run ID for an autosomal migration.

`00_metadata/chromosome_policy.json` records the frozen policy/coordinates; `00_metadata/masks/RES.scope.json` records actual per-matrix scope and registration. Inherited weights retain their upstream balancing scope. Common downstream coordinates do not imply that upstream weight calculations used identical chromosomes. Interpret processing differences explicitly.
