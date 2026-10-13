# General-purpose upgrade: microc2compare 0.4.0

The new modules operate on the current project's sample manifest, registered
reference assets, native matrices and frozen run metadata. A completed mouse
comparison, a prior loop-export bundle and particular observation names/counts
are not prerequisites. Arbitrary reference identifiers are supported; their
chromosome lengths and view scope must be supplied. Assembly-matched annotations
are required only when annotation is enabled.
For browser export, the chromosome-size file must include every native matrix
chromosome, including scaffolds excluded from analysis views. Native coordinate
profiles require these lengths even for UCSC-only output.
The example `mm39` default is a convenience, not an assembly restriction.

## Start a new project

```bash
microc2compare --init /owned/human/config/project.conf --reference hg38 --preset practical_2_4kb
# Edit the manifests/reference paths and experimental identity before running.
# Optional resource settings, sourced from the project configuration:
# source /installed/microc2compare/config/resources_72cpu_400gb.example.conf
microc2compare -c /owned/human/config/project.conf --prepare-inputs
microc2compare -c /owned/human/config/project.conf --register-resources
microc2compare -c /owned/human/config/project.conf --input-capabilities
microc2compare -c /owned/human/config/project.conf --audit-resolutions
microc2compare -c /owned/human/config/project.conf --dry-run
microc2compare -c /owned/human/config/project.conf
```

Preparation actions return new configuration variants; use the configuration
they print. They are explicit actions, separate from scientific execution.
Imported counts and stored weights remain read-only. Existing configurations
retain their defaults; new `--init` projects enable automatic UCSC export and
set an eight-observation admission ceiling. CPU/RAM budgets can admit fewer.
The 72-CPU/400000-MB profile applies only when the actual host/allocation permits
it. It reserves at most eight concurrent observations, up to eight workers per
caller and eight workers in run-level aggregation/quantification. A ceiling is
not a promise of eight useful workers for every grid or stage.

## CPU execution and interruption recovery

The Bash scheduler fills a freed observation slot immediately, rather than
waiting for a whole wave. Managed task process groups are stopped on failure or
interruption. QC has independent per-grid workers; Mustache and Chromosight use
bounded internal workers; TAD calling retains its existing per-observation
parallel implementation. Raw count extraction and APA/ADA/saddle measurements
use independent processes that open their own read-only matrix handles.
Native BLAS/OpenMP pools are held to one thread to avoid nested oversubscription.

Stage reservations, observation ceilings, effective CPU affinity/cgroups and
available RAM bound admission. Array-based memory estimates are conservative
planning estimates, not kernel-enforced RSS limits. Inspect actual resource
records before raising concurrency. `LIMIT_MATRIX_READERS=true` adds a
conservative total reader limit using reserved worker counts; lower stage thread
reservations when `MAX_PARALLEL_IO` is smaller. Disk throughput may saturate
before CPU occupancy does.

Each aggregate observation/module/grid publishes atomic products and a
checksum-verified unit checkpoint. Interrupted units are recomputed; valid
units can be reused even if worker counts change. The parent merges committed
tables in manifest order and only then writes condition summaries. Expected
curves are generated from fresh inputs and cached with matrix identity,
resolution, native mask, labels, weights, views, diagonal policy, implementation
and runtime/environment identities. A previous cohort's cache is unnecessary.
`EXPECTED_CACHE=false` or `AGGREGATE_CHECKPOINTS=false` disables reuse.

Completion validation also reuses prerequisite checks and source hashes within
one audit. This ephemeral cache expires when the public validation call returns;
the next call checks inputs, outputs, dependencies and code afresh. It removes
repeated traversal of the same dependency graph without accepting stale records.

## Fine discovery and approximately 4-kb separation

`hsieh_8grid` requests 400, 600, 800, 1000, 2000, 4000, 10000 and 20000-bp grids;
`short_range_4kb` requests the first six and a shorter maximum separation.
Both are adapted CPU strategies, with Mustache and an independent Chromosight
pattern-discovery stage. Matching grid sizes does not reproduce the historical
Hsieh analysis. See [the caller profile and exact-pairs documentation](cpu_loop_profiles.md).

An element 4 kb from a promoter needs sufficiently fine bins and a compatible
caller background. A 4-kb *bin* does not meet that separation goal. The tested
Mustache settings have a best-case six-bin occupancy floor: 2.4 kb at 400 bp,
3.6 kb at 600 bp, 12 kb at 2 kb, and 24 kb at 4 kb. Masks, sparse data and other
filters can make the observed limit larger. Synthetic recovery of a 4-kb peak
tests capability, not sensitivity or specificity in a biological study.

Coarse maps cannot supply authentic finer grids. `--prepare-resolutions`
creates exact integer coarsenings only. `--prepare-from-pairs` requires complete,
deduplicated retained exact-position pairs with explicit source SHA256,
coordinate system, filters and balancing recipe. It creates new grids, checks
raw contact conservation and publishes a registry only after convergence.
The scientific pipeline never performs this preparation implicitly.

Sub-kilobase measurements require `ALLOW_SUBKILOBASE_ANALYSIS=true`. Loop
measurement/discovery at 1 kb or below also requires a project-specific reviewed
resolution audit, including contact retention and reproducibility. The audit
action emits `eligibility.template.json` with **no approval**; record the actual
review and point `RESOLUTION_ELIGIBILITY_REGISTRY` to that completed registry.
Source-call membership retains caller and resolution identity; support counts
distinct observations, not the number of software tools detecting a pixel.
Caller enrichment or exploratory correlation P values are not differential
biological significance. Formal comparisons retain independent-feature,
replication, design and calibration gates.

## Stripes

`RUN_STRIPES=true` adds a separate optional StripePy 1.3.0 CPU stage and
descriptive cross-observation support. It requires a new `STRIPE_ENV_PREFIX`
using [envs/stripes.yml](../envs/stripes.yml); the core environment is unchanged.
The tested adapter streams raw cis contact belts into task-owned temporary maps
split at common usable-bin gaps. It detects only within eligible segments;
cross-gap features are unavailable. Original chromosome coordinates, seed,
orientation, native descriptors, segment mappings and exclusion reasons are
retained. No-call output is a validated empty result with coverage information.
Native relative change is an absolute percentage measured on StripePy's
transformed signal; retained calls also require inner signal above outer signal.
It is not the raw loop O/E ratio or a condition effect.
See [stripe detection](34-stripe-detection.md).

## Unified browser tracks

`RUN_BROWSER_EXPORT=true` writes one dated `run/browser_tracks/` bundle after
enabled analyses and before the HTML report. It includes individual observations
and condition summaries, native loop calls and anchors, available TADs/boundaries,
exact ADA domains/scalars, all analysed insulation windows, compartment
components actually saved, selected oriented scores/A/B assignments and stripes.
Missing components/modules are recorded explicitly. Annotation is optional.

Continuous condition means require compatible coordinates and the necessary
orientation; missingness and contributions are retained. Discrete tracks use
documented detection/consensus support rules. They never average coordinates or
claim significant specificity. ADA conditions retain matching stable domain IDs.
Compartment means use selected phased scores only; raw E1/E2/E3 components are
exported per observation without condition averaging. Their saved scaling is
preserved, so mean amplitude is not a quantitative compartment-strength estimate.
The README, track inventory, source identities and checksums describe loading.
UCSC files contain the appropriate custom-track headers.

```bash
# Export saved results independently, without rerunning scientific analyses:
microc2compare -c /owned/project.conf --export-browser
microc2compare -c /owned/project.conf --export-browser --out /owned/new-browser-bundle
```

HiGlass export is optional: `EXPORT_HIGLASS=true`. Continuous signals use
bigWig; domains/boundaries use interval databases; loop/stripe rectangles use
two-dimensional databases. Loops are not coerced into bigWig. Set a separate
`BROWSER_ENV_PREFIX` from [envs/browser.yml](../envs/browser.yml) and use explicit
bootstrap to install/check converters. Without available converters the inventory
records missing conversions rather than claiming files were generated.
Ordered chromosome profiles distinguish matrix offsets even within one assembly.
Conversion, server ingestion and actual rendering have separate statuses; the
workflow does not publish or ingest files. See [browser export](browser_export.md).

## Observation identity and migration

`subject_id` remains the donor/independent organism identifier. Optional
`experimental_unit_id`/`culture_id` supports independently grown cultures from
one donor for descriptive analysis. Repeated donor/condition rows need distinct
explicit culture units. Two cultures are not two donors. The existing R model
does not handle nested cultures; inferential nested designs are rejected rather
than treated as extra independent replication. Reserved observation names include
`mask`, `labels`, `pool__*` and `condition__*`.

Install 0.4.0 alongside the accepted older release. Use new run IDs for scientific
runs; do not replace code in a running task or manufacture new completion records
for old results. Exporting validated saved tables is independent from rerunning
their analyses. Existing completed biological results and shared server aliases
are not altered by this source upgrade.

## Validation boundaries

Local checks include fresh projects with unequal condition sizes and arbitrary
assemblies, source checksum preservation, serial/parallel numerical agreement,
checkpoint recovery, Mustache geometry/band access and actual Chromosight APIs and an actual one/two-process
StripePy call in an isolated runtime. Linux CI definitions exercise the Bash
scheduler, optional prefixes and converter paths. Remote CI, shared-server
performance, real-data short-range/stripe sensitivity, browser rendering and
formal biological calibration remain separate acceptance steps. Source tests
do not certify all of those steps.
