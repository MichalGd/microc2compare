# Generic CPU loop profiles and input capability reports

`loop_profiles.profile_defaults(name)` returns fresh configuration values for
`hsieh_8grid` and `short_range_4kb`. It supplies no assembly, sample identifiers,
replicate count or prior pilot data. The eight-grid preset requests 400, 600,
800, 1000, 2000, 4000, 10000 and 20000 bp. The short-range preset requests the
first six grids with a 500 kb maximum separation. Both retain scientific
admission for <=1000 bp grids; neither grants coverage or reproducibility
approval. The 120000 bp APA flank divides every requested grid.

The eight grids, Mustache p threshold 0.1, sparsity 0.88 and two octaves are
documented in the primary [Hsieh et al. 2022 paper](https://www.nature.com/articles/s41588-022-01223-8)
and [supplementary methods, printed pages 8–9](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41588-022-01223-8/MediaObjects/41588_2022_1223_MOESM1_ESM.pdf).
These presets adapt that study: the tested software versions are Mustache
1.3.3 and Chromosight 1.6.3, normalization uses registered upstream weights and
common-view O/E, tiled discovery replaces the historical CLI, and merging uses
the project's explicit complete-linkage correspondence. Historical caller
distances, thresholds and depth are not claimed as a 4 kb detection recipe.
The original `hsieh_inspired` profile retains its source-registry gate.

`resolutions.input_capabilities(c, samples, requested=None)` is a read-only
metadata report. Each sample/grid is available, invalid, explicitly coarsenable
or missing. Retained-pairs availability is reported separately from authenticity
or adequacy. A 1000 bp matrix cannot supply authentic 400, 600 or 800 bp pixels.
All eight grids can be derived from an authentic 200 bp base, because it divides
every target grid, or directly from retained exact-position deduplicated pairs.
Coarsening needs a new integer-count map and newly estimated weights with an
explicit recipe. Coarse weights are never copied as fine-grid weights.

`mustache_geometry(c, resolution)` describes the tested caller's necessary ideal
occupancy bound. The first assigned scale is `sigma_zero * 2**(2/10)`. At sigma
1.6 its inner radius is two bins. The original nonzero mask requires diagonal
>=4: at a five-bin offset only 19/25 inner pixels are geometrically possible,
failing sparsity 0.88. At six bins 22/25 pass, and the outer 60/81 window passes
the 0.6 threshold. Thus the ideal floors are 2400/3600/4800/6000/12000/24000/
60000/120000 bp across the eight grids. At 800 bp a nominal five-bin 4 kb offset
is ineligible with these parameters. Zeros, masks, boundaries, larger scales,
enrichment and significance can increase the observed floor. This derivation
uses the [tested Mustache source](https://github.com/ay-lab/mustache/blob/master/mustache/mustache.py);
the [Mustache paper](https://doi.org/10.1186/s13059-020-02167-0) is a Genome Biology
article, not DOI 10.1038/s41467-020-18507-9.

`mustache_resource_plan` limits workers by CPU and the configured memory budget,
reserving 20% overhead and charging known global sparse arrays separately.
Mustache tiles have side `max(2 * floor(max_distance/resolution), 2000)`; the
estimate reserves twenty float64 arrays per worker. At 400 bp/2 Mb one worker
needs an estimated 16 GB of dense arrays. A configured 80 GB reservation initially
admits four workers, dropping to three after charging nonempty shared arrays.
Actual per-view plans are returned in `mustache_calls(...).attrs['resource_plans']`.
The upper band is read in bounded row blocks rather than as a chromosome-wide
sparse square. These estimates are conservative admission heuristics, not
measured RSS caps; representative project data still needs resource validation.

`fine_loops.chromosight_discovery(c, scope)` independently emits
`03_features/loops/<scope>.<grid>.chromosight.tsv` and matching JSON provenance.
It reads the registered matrix, applies inherited weights and common-view O/E,
and uses the real Chromosight normalized correlation and native focus picker.
Rectangular tiles bound memory; the first-anchor core row owns each candidate.
Kernel halos overlap, but focus components are grouped within each tile, so
very broad components spanning tiles can differ from chromosome-wide CLI
grouping. Missing-pixel tolerance is configured, observed zeros must be below
10%, and candidate centres obey the requested genomic distances. Pixels around
the centre provide kernel context beyond that centre-distance interval.
Calls preserve measured native coordinates; no caller snaps a peak to a desired
4 kb offset. The [Chromosight paper](https://www.nature.com/articles/s41467-020-19562-7)
describes the underlying pattern detector.

Discovery is CPU only and never invokes Chromosight's automatic matrix
normalization. Per-window caller p-values are exploratory; caller q-values and
biological significance remain unavailable. Independent discovery from the
same contact maps is not independent biological validation. Fixed-candidate
scoring remains a separate function. Its threads use separate lazy read-only
Cooler objects, preserve candidate order, and reuse expected values only from a
hash-verified completed QC stage for the current run. Score tables retain missing
values for ineligible, masked, flat and edge windows.

The dedicated tests use real Cooler reads for Mustache bands and the installed
Chromosight kernel/correlation, multiple sample counts and assemblies, missing
fine grids, punctuated sample IDs, and unchanged source SHA256 hashes. A
synthetic 4 kb peak is recovered within one 400 bp bin; this checks software
behavior and is not a depth or sensitivity claim for a biological dataset.

## Explicit exact-position pairs preparation

`pairs_preparation.prepare_from_pairs(c, out, resolutions=None)` runs only when
explicitly requested. It requires a new output directory and complete biological
sample/technical-library metadata. Each sample supplies `pairs_path` and this
record inside its provenance JSON (replace the hash and declared values with
the actual upstream evidence):

```json
{
  "reference_id": "your_assembly",
  "pairs_provenance": {
    "reference_id": "your_assembly",
    "source_sha256": "SHA256_OF_RETAINED_PAIRS_FILE",
    "complete": true,
    "deduplicated": true,
    "dedup_scope": "consolidated_biological_unit",
    "coordinate_system": "1-based",
    "upstream_filters": {
      "mapq": 30,
      "chromosome_policy": "all_registered_chromosomes",
      "retained_pair_types": ["UU"],
      "minimum_cis_distance_bp": 1000
    },
    "balancing_recipe": {
      "cis_only": false,
      "ignore_diags": 2,
      "mad_max": 5,
      "min_nnz": 10,
      "min_count": 0,
      "tol": 0.00001,
      "max_iters": 1000,
      "rescale_marginals": true
    }
  }
}
```

The illustrated balancing recipe is an explicit declaration, not an assurance
that these thresholds suit a particular project. `dedup_scope` can also be
`joint_library`. Deduplication is declared upstream and is not recomputed here.
The pairs declarations must include their actual mapping-quality threshold and
chromosome policy. Those declarations, the retained-pairs distance cutoff and
the validated technical membership supply the NEW maps' active scientific
provenance. An older matrix's filter cutoff and balancing recipes are preserved
under `source_provenance`; they never replace the NEW maps' cutoff or recipes.
The pairs file needs a unique `#columns:` header containing `chrom1`, `pos1`,
`chrom2`, `pos2` and `pair_type`. A `#genome_assembly:` value, when present, must
match the registered reference. Plain text and gzip/BGZF text are supported.
The reference registry must bind complete chromosome sizes by SHA256.

Positions are validated as exact integers and converted using the explicitly
declared zero- or one-based convention. Chromosomes and bounds must match the
registered sizes. Declared pair types are retained; unmapped `!`/`*` anchors are
excluded and counted. A contact below the declared upstream cis-distance filter
is an error, since it contradicts the record. Finer grids cannot restore contacts
already removed upstream. Optional `records_total` verifies completeness against
the parsed record count.

Each requested grid is counted directly from the same retained coordinates,
with streaming aggregation and Cooler's external merge. This scans the source
once per requested grid and may use substantial temporary disk space. Every raw
integer map must conserve the retained contact count, and its new balancing must
converge under the supplied recipe. Missing weights produced by that recipe
remain missing. No source map or weight is copied or modified. Before the usable
registry is published, every entry is checked for exact registered bins,
integer-count conservation, valid converged weights, full scientific provenance
and technical membership. A failure
leaves a diagnostic and unfinished NEW maps for inspection without publishing a
usable registry.

After all grids succeed, the returned dict names `matrix_registry`, `samplesheet`
and `outputs`. Set the project to use those new registrations, then run preflight
and the project-specific fine-grid retention/coverage/reproducibility review.
Preparation publishes `scientific_admission: not_performed` and never approves
the resolution eligibility registry.

The Windows development runtime lacks Mustache's optional `hicstraw` import and rejects a mandatory dual-caller probe explicitly. Hosted Linux acceptance at commit `8605a236ab17` passed both actual API probes, CPU/pairs tests and full/selected fine workflows. Later generic source corrections have additional focused checks; the current branch CI repeats the required APIs and complete fine acceptance. These compatibility fixtures do not calibrate biological loop sensitivity.
