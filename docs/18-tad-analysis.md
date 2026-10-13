# TAD analysis and descriptive changes

A topologically associating domain (TAD) is a region with preferential internal contacts. This module uses **insulation-defined domains**: intervals between consecutive strong insulation boundaries within a registered chromosome/view. It does not reconstruct nested/hierarchical domains. Window size, resolution, coverage and boundary threshold influence which domains are detectable.

## Inputs and running

Use compatible biological-unit matrices at `INSULATION_RESOLUTION`, upstream stored weights, a samples table, technical membership, registered chromosome sizes/views and optional blacklist. A phasing track and the R environment are not needed. `comparisons.tsv` declares case/reference conditions; for a project without comparisons use a header-only table. Calling domains does not create biological replication.

These are **Linux Bash** commands. Set paths before using them:

```bash
m2c_src="/opt/conda_envs/microc2compare/releases/agreement-20261008/microc2compare-0.2.3"
m2c_conf="$HOME/microc2compare-work/projects/pilot/config/project.conf"
bash "$m2c_src/scripts/microc2compare.sh" -c "$m2c_conf" --analysis tads --dry-run
bash "$m2c_src/scripts/microc2compare.sh" -c "$m2c_conf" --analysis tads
```

The selection enables insulation and necessary input registration, then TAD calls/events and a report. It does not enable loops, compartments, HiCRep or R comparisons. For an existing run use the same command with `--resume`. In the configured full pipeline, `RUN_TADS=true` is active when `RUN_INSULATION=true`.

## Parameters

Put overrides in your project `.conf` after any resource preset. They are snapshotted in the run. Defaults are initial descriptive policies to evaluate in your biological pilot, not universal biological cutoffs.

| Parameter | Default | Meaning and effect |
|---|---|---|
| `INSULATION_RESOLUTION` | `10000` bp | Existing grid; determines coordinate precision |
| `INSULATION_WINDOWS_BP` | `50000,100000,250000` bp | Scales at which crossing contacts are measured; primary window must be included |
| `INSULATION_PRIMARY_WINDOW_BP` | `100000` bp | Primary report/insulation comparison scale; other TAD windows are retained separately |
| `INSULATION_BOUNDARY_THRESHOLD` | `Li` | Existing cooltools automatic boundary-strength threshold; `Otsu` or a finite fixed number is supported |
| `INSULATION_MIN_VALID_PIXEL_FRACTION` | `0.8` | Minimum eligible/full-diamond coverage; also used for domain/evaluable-region score coverage |
| `INSULATION_IGNORE_DIAGS` | `2` bins | Ignore near-diagonal contacts in the calculation |
| `TAD_MIN_LENGTH_BP`, `TAD_MAX_LENGTH_BP` | `30000`, `3000000` bp | Accepted boundary-to-boundary lengths |
| `TAD_BOUNDARY_TOLERANCE_BP` | `10000` bp | Allowed complete-span coordinate spread for condition consensus and matching |
| `TAD_MIN_REPLICATE_FRACTION` | `0.67` | Required domain/evaluable-region support; round `fraction * number_of_samples` upward. With three samples, this default requires all three; use a reviewed alternative if two-of-three is intended |
| `TAD_MIN_RECIPROCAL_OVERLAP` | `0.5` | Required one-to-one shared fraction of both intervals; smaller-over-larger overlap admits candidate split/merge connections |
| `TAD_SIZE_CHANGE_BP` | `20000` bp | Absolute length change must exceed this threshold |
| `TAD_SIZE_CHANGE_FRACTION` | `0.1` | Fractional change relative to reference length must also exceed this threshold |
| `TAD_STRENGTH_CHANGE` | `0.25` | Absolute mean endpoint-boundary prominence change used to annotate otherwise stable matches |
| `TAD_SENSITIVITY_THRESHOLDS` | `0,0.2,0.5` | Additional fixed thresholds evaluated on the same measured scores; no extra matrix copies |

Li/Otsu select strong boundaries from the prominence distribution; numeric thresholds are in the tool's boundary-strength units, not a P value or a fold change. The numeric rule includes values equal to the threshold. Alternative automatic methods or windows require a separate insulation variant; the within-run sensitivity product varies fixed thresholds on existing window scores.

## Calls, matching and uncertainty

Each biological sample is called separately. Domains cannot bridge excluded bins. Condition-domain clustering bounds the complete spread of both endpoints by the tolerance and allows each biological observation one vote per cluster. A real observed interval is the deterministic representative, rather than an averaged coordinate unsupported by any sample. Unsupported and conflicting overlapping condition intervals are retained in a disagreement table.

Within a comparison, construct correspondences only within the same chromosome/view. Coherent adjacent many-to-one and one-to-many matches become merges/splits; a one-to-one match must pass reciprocal overlap. Each connected event group is classified once. A merge's absorbed domains are not also counted as independent disappearances.

| Label | Meaning |
|---|---|
| `merge` | Adjacent reference domains match one case domain with compatible outer boundaries |
| `split` | One reference domain matches adjacent case domains with compatible outer boundaries |
| `growing`, `shrinking` | A one-to-one domain changes length beyond both absolute and fractional thresholds |
| `shifted` | Boundary displacement exceeds tolerance without an accepted large length change |
| `appearing`, `disappearing` | Detectable domain gain/loss without a merge/split explanation, with adequate data in both conditions |
| `strength_change` | Otherwise stable geometry has endpoint-prominence change above threshold |
| `stable` | No change beyond the declared geometry/strength thresholds |
| `complex` | Many-to-many or insufficiently coherent correspondence |
| `uncertain` | Insufficient common coverage, view-edge evidence, or overlapping unsupported/conflicting replicate calls |

Direction is **case relative to reference**. Growing/shrinking refers to genome length in bp, not stronger/weaker contacts. The coverage check requires both outer boundary bins to exist and be usable, adequate internal score coverage and an uninterrupted common mask in the required fraction of biological samples. The current event check assesses measurement eligibility; it does not estimate a probability of biological absence. Genuine gains/losses near weak boundaries still require sensitivity and locus review.

Endpoint prominence and mean log2 insulation are recorded separately. Prominence changes in multi-domain groups summarize participating endpoints, which may differ; they are not a test at one fixed boundary. Intra-domain O/E contact strength is separately available from ADA on the fixed consensus domain universe. Boundary strength, insulation score, O/E and domain length are different measurements.

## Exact outputs

Paths below are relative to the run root. `SAMPLE` is a safe biological sample ID; `WINDOW` is bp. Condition filenames use a deterministic `condition_HASH` identifier because user condition labels can contain spaces; `parameters.json` maps IDs back to labels.

| Path | Content |
|---|---|
| `03_features/tads/samples/SAMPLE.WINDOW.tsv` | Per-replicate domain coordinates, endpoint prominence, mean insulation, length, sample ID and descriptive status |
| `03_features/tads/samples/SAMPLE.WINDOW.bed` | Browser intervals: chromosome, start, end, feature ID |
| `03_features/tads/samples/SAMPLE.WINDOW.exclusions.tsv` | Candidate intervals rejected by length, masks, boundary eligibility or score coverage; reason is explicit |
| `03_features/tads/conditions/condition_HASH.WINDOW.tsv` | Condition-supported intervals with support count/fraction and contributing sample IDs |
| `03_features/tads/conditions/condition_HASH.WINDOW.unsupported.tsv` | Insufficient-support or conflicting intervals; inspect replicate disagreement |
| `03_features/tads/conditions/condition_HASH.WINDOW.bed` | Condition-supported browser intervals |
| `06_comparisons/ID/tad_events.WINDOW.tsv` | Event group, participating reference/case domain IDs, coordinates, sizes, shifts, support, coverage and reason |
| `06_comparisons/ID/tad_events.WINDOW.bed` | Affected event-group intervals for browsers |
| `06_comparisons/ID/tad_events.WINDOW.summary.tsv` | Counts of primary event groups by type |
| `03_features/tads/sensitivity/ID.WINDOW.THRESHOLD.events.tsv` | Complete event classifications at each additional fixed threshold |
| `03_features/tads/sensitivity.tsv` | Counts by comparison, window, threshold and event; zero-event settings are recorded |
| `03_features/tads/parameters.json` | Method, coordinate convention, condition-ID mapping and settings |
| `08_reports/figures/tad_events.png`, `tad_lengths.png` | Available primary-window event counts and matched size-change distribution |
| `08_reports/figures/locus.SAMPLE.png` | Bounded contact-map examples with per-replicate domain boundaries; long domains may extend beyond the displayed region |

All BED coordinates are zero-based, end-exclusive. TSV headers define columns. For example, an event with `reference_length_bp=200000`, `case_length_bp=300000`, `length_change_bp=100000`, and `length_ratio=1.5` is an expansion if matching/coverage and configured size thresholds pass. P/FDR columns remain NA. A domain disappearance has missing case coordinates/length rather than an invented length of zero.

Header-only TSVs and comment-only empty BEDs are valid zero-call products. Empty calls do not mean there are no biological domains; review scale, threshold, masks and coverage. A computation failure is never converted to an empty success. Browser loading may require selecting the correct assembly; event tracks are annotations, not significance tracks.

## Parameter variants and validation

Use [the variant helper and commands](17-analysis-modules.md#reports-variants-and-errors) to preserve earlier analyses. Compare event tables and bounded loci across windows/thresholds, not only totals. Fixed-threshold sensitivity can change domain identities as well as counts. Short boundary movements at the grid/tolerance scale should not be interpreted as precise biological relocation.

The tests cover known interval events, direction, masks, thresholds and disagreement. The scientific integration test additionally exercises actual cooltools insulation and the TAD adapter when that runtime is available. The new candidate needs target-server acceptance and a biological pilot. These classifications have **no event-level biological P values/FDR**. Calibrated fixed-grid differential insulation supplies related evidence but does not transfer bin-level FDR to a domain or event. A formal TAD-event model is a separate future method.

Method context: [cooltools insulation](https://cooltools.readthedocs.io/en/latest/notebooks/insulation_and_boundaries.html) and [TADCompare's primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7076128/). This module uses its own insulation/interval rules; it does not run TADCompare or reproduce its statistical model.

## Execution and performance

Version 0.3.3 reuses view-local indexes across thresholds and coverage/disagreement queries. Independent sample domain calls use up to `THREADS_TADS` processes (default 2); `MEMORY_MB_TADS` defaults to 8000 decimal MB. These defaults are internal, so existing configs need no new fields. Optional larger reservations stay within total/effective CPU and estimated memory bounds. The parent preserves deterministic matching and output ordering.

`03_features/tads/performance.json` records actual workers and timings. See [performance and parallelism](23-performance-and-parallelism.md) for tuning, numerical agreement tests, and acceptance that replays an earlier completed pilot into a new directory without modifying it.
