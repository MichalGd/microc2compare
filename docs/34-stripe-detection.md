# Optional CPU stripe detection

Stripe detection is optional (`RUN_STRIPES=false` by default). The adapter uses
[StripePy 1.3.0](https://stripepy.readthedocs.io/en/stable/) on raw counts and
records descriptive geometric and signal measurements. Install the separate
[stripe environment](../envs/stripes.yml) and configure its prefix through the
software profile; its NumPy 2 dependencies are separate from the accepted loop
environment. The supported version is checked before discovery.

The recipe pins Python 3.11, NumPy 2.2.6 and SciPy 1.17.1. The
[conda-forge SciPy package](https://anaconda.org/conda-forge/scipy) provides
Python 3.11 builds compatible with this NumPy version. SciPy 1.17.1 was also
used in the local real-caller tests, with NumPy 2.4.6; the exact Linux recipe
is checked separately by the stripe CI job.

```bash
microc2compare -c /your/project/project.conf --analysis stripes --dry-run
microc2compare -c /your/project/project.conf --analysis stripes
```

Set `STRIPE_RESOLUTIONS="4000"` or an explicit comma-separated list such as
`"2000,4000"`. Every observation needs the requested grids in its original map or
registered matrix overrides. Assembly, chromosomes, views, masks, conditions,
replication and resource budgets come from the project; no completed cohort is
needed. The adapter does not recalculate imported weights or need an expected
curve. `STRIPE_NORMALIZATION="none"` denotes raw-count discovery.
Inherited balancing metadata and finite positive weights remain required to
define common eligible bins; the caller receives raw pixels without those
weights. An excluded inherited-weight bin remains a gap in a stripe-only run.

The current `STRIPE_MASK_POLICY="common_segments"` splits each registered view
at every excluded bin. A task-owned temporary Cooler contains only raw cis
contacts within these contiguous usable segments and the requested genomic
belt. The adapter validates integer count conservation, removes this temporary
map after calling, and keeps the source maps unchanged. Segments shorter than
three bins are recorded as too short. The caller receives `--min-chrom-size 1`,
so its usual 2 Mb chromosome cutoff cannot silently discard the segments. This
conservative analysis excludes stripes crossing gaps: no retained calls does
not establish biological absence or complete genome coverage.

`THREADS_STRIPES` sets caller CPU processes. `STRIPE_GENOMIC_BELT_BP` and
`STRIPE_MAX_WIDTH_BP` are physical lengths. The pinned caller supplies its other
published defaults. Local CPU time and memory must be measured; the published
10 kb runtime does not predict 2/4 kb runtime or accuracy. The
[StripePy paper](https://doi.org/10.1093/bioinformatics/btaf351) benchmarks
synthetic maps at 5/10/25/50 kb and real maps at 10 kb. Finer grids require depth,
replicate and visual validation.

Each `03_features/stripes/SAMPLE.RES` prefix contains:

- `.hdf5` and `.caller.tsv`, retaining the native caller result when eligible
  segments exist. Their coordinate axis is segment-local, not original genomic
  coordinates. Use `.segments.tsv` to recover original offsets.
- `.tsv` and `.bedpe`, canonical original-genome intervals, with the narrow
  anchor first and extended interaction body second. Native inclusive HDF5 bin
  bounds become zero-based half-open BED coordinates by adding one bin to each
  upper bound; the unmodified bounds are also retained in the TSV. Consequently
  the official caller's BEDPE view and the canonical BEDPE have separately
  documented endpoint conventions.
- `.segments.tsv`, `.exclusions.tsv` and `.json`, describing selected scope,
  unusable segments, rejected predictions, counts, normalization, software and
  input identities, commands and coordinate conversion.

StripePy 1.3.0's `view`/descriptor reader raises an assertion for a triangle with
zero predictions. The adapter reads its geometry first and skips only genuinely
empty triangles before fetching their biological descriptors. `.caller.tsv`
therefore exports native inclusive-bin records directly from the HDF5 instead
of invoking the failing `view` command. Other caller/read errors remain failures.

The body includes the diagonal seed bin, so it may overlap the anchor. `LT`
corresponds to a downstream genomic excursion and `UT` to an upstream excursion
under the caller's native matrix convention. These directions are independent
of gene strand. Matrix rectangles and separate anchor/body BED tracks convey
stripe geometry; ordinary focal-loop arcs do not.

In StripePy 1.3.0, `relative_change` is the native percentage
`abs(inner_mean-outer_mean)/outer_mean*100`, computed on its transformed signal.
The adapter requires positive contrast (`inner_mean > outer_mean`), finite
contrast and a positive outer mean. `STRIPE_MIN_RELATIVE_CHANGE` sets the
percentage cutoff; zero retains every otherwise eligible enriched prediction.
These values are not contact-count fold changes or biological P values/FDR.

`summary.tsv` records each observation's call count and its condition's actual
observation denominator, together with a separate distinct-subject count.
Independent culture labels within a shared donor do not create additional
donors. `matches.tsv` retains pairwise, potentially many-to-many
anchor-overlap matches with the same genomic excursion direction at the same
resolution, together with body-overlap measurements. `summary.json` records
dynamic condition membership. Caller agreement, multiple resolutions and
technical inputs do not create biological replicates. Statistical differential
stripe inference remains deferred.

The report preserves no-call observations and unavailable current results in
`08_reports/stripe_observations.tsv`, and shows eligible-segment coverage as a
percentage of selected registered-view length. Condition counts, direct stripe
support, distinct-subject counts, native relative-change ranges and a call-count
chart remain descriptive. A singleton's 1/1 support does not supply replicate
confirmation. The report reads only products validated for the current stage
signatures; saved stale tables are not treated as current evidence.

The adapter tests separate mocked interface tests from an optional real
StripePy CLI tests with one and two CPU processes on a new masked synthetic map.
The real tests run only in the
pinned stripe environment; passing it establishes software compatibility,
coordinate/mask handling and unchanged inputs, not calibration of biological
stripe calls. Linux release acceptance must also exercise interruption, CPU
admission, resumability and report/browser consumption.
