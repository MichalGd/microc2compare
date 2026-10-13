# Export browser tracks for a new project

Browser export reads the current run's frozen configuration, sample metadata, registered chromosome grids and saved result tables. It does not depend on an earlier project or earlier browser export. Reference identity, chromosome names/order, conditions and observation counts come from that run. Any assembly with a valid registered reference is supported; chromosome labels are never renamed or lifted over.

The registered chromosome-size file must cover every chromosome on each native matrix axis, including scaffolds excluded from the analysis views. The exporter builds native coordinate profiles even for UCSC-only output, so a reference containing only selected analysis chromosomes is insufficient when the matrices contain other chromosomes.

Enable the final workflow stage in the project configuration:

```bash
RUN_BROWSER_EXPORT="true"
EXPORT_UCSC="true"
EXPORT_HIGLASS="false"
BROWSER_OUTPUT_DIR=""
HIGLASS_CHROMSIZES=""
```

An empty output directory setting creates `RUN_DIR/browser_tracks/<UTC timestamp>/`. An explicit `BROWSER_OUTPUT_DIR` must name a new directory. Existing bundles are never overwritten. These export options do not select scientific modules: enable the analyses needed by the project separately. Missing optional modules are recorded in `missing_features.tsv`; an incompletely saved measurement module fails export instead of manufacturing condition means.

Export saved results independently, using the installed repository's script:

```bash
python scripts/export_browser.py --run /path/to/current/run
python scripts/export_browser.py --config-json /path/to/resolved/config.json --out /path/to/new/bundle
python scripts/export_browser.py --run /path/to/current/run --higlass
```

The main CLI also accepts `microc2compare -c PROJECT.conf --export-browser --out /path/to/new/bundle` for saved results; the action performs no scientific stage reruns.

Use `--reference-registry` when registered reference assets have moved; their reference ID and chromosome-size checksum must still match. Preserve `00_metadata/inputs.json` when copying a run: the exporter checks the original registry against its frozen fingerprint when available, and otherwise requires the relocated chromosome-size file to match a frozen input fingerprint. A newly declared checksum in a replacement registry is insufficient. `--higlass-chromsizes` selects an explicit ordered profile with matching registered lengths that covers all exported chromosomes. `--no-ucsc --higlass` retains raw tables and geometry while omitting UCSC custom headers. Exporting saved tables requires no contact matrices or scientific caller reruns.

Each bundle has one `README.md`, `inventory.tsv`, `samples.tsv`, `source_files.tsv`, `artifact_checksums.tsv` and `export_audit.json`. Products use `ANALYSIS/individual/SAMPLE/` or `ANALYSIS/conditions/CONDITION/`, with shared products under `ANALYSIS/shared/`. UCSC files are in `ucsc/`; converted HiGlass files and domain/stripe rectangle geometry are in `higlass/`. Sidecars preserve stable feature IDs, original scalar values, support counts and the source-table fields. `reference/coordinate_profiles.json` and `profile_assignments.tsv` identify ordered chromosome names and lengths for the registered analysis and native matrix axes. Profile paths are relative to the bundle.

Upload nonempty `*.custom.txt` files to UCSC Custom Tracks with the matching genome assembly. The inventory distinguishes custom tracks, raw BED, bedGraph and BEDPE, rather than treating a TSV table as a browser format. All coordinates are zero-based half-open.

| Saved result | Exported meaning |
| --- | --- |
| Loops | Native per-grid calls retain their original endpoints. Catalogue representatives are separate; sample and condition support count each original observation once. Annotation is optional. UCSC interact and true endpoint BEDPE are separate formats. |
| Compartments | Selected oriented scores, raw E1/E2/E3, component identity and sign agreement remain distinct. Condition means use selected phased scores only; raw components are exported per observation without condition averaging. |
| Insulation/boundaries | Every saved/configured window is exported with exact bins, finite-value coverage, individual calls and saved consensus support. |
| TADs | Exact sample and condition intervals, unsupported/excluded records and saved supporting-domain means are retained. Each domain is also an exact square on both matrix axes. |
| ADA | Exact primary-window domain IDs and mean intra-domain O/E, sums and eligible-pixel counts are retained. A domain scalar spans that interval; it is not a per-bin signal or proof of resized pileup membership. Overlapping domains retain BED/rectangle/TSV records; overlapping scalar bedGraph is omitted. |
| Stripes | Individual and condition tracks retain every exact native anchor/body rectangle, including diagonal seed overlap, direction, triangle and descriptors. Condition tracks annotate direct support from distinct observations; no merged consensus geometry or differential significance is inferred. |

New compartment, insulation and ADA condition means give each original observation equal weight and require complete finite coverage. Denominators come from metadata, including unequal condition sizes. Available-case summaries and finite counts are separate. Saved TAD means retain their original supporting-domain semantics. Detection support is descriptive and does not establish biological absence.

Selected compartment components can differ between views and observations. The exporter preserves their saved scaling without further alignment, standardization or rescaling. Raw component signs remain arbitrary; a descriptive condition mean is neither a newly pooled eigenvector nor a quantitative compartment-strength estimate.

For each condition and stripe grid, `direct_support_<grid>.tsv` retains original stripe IDs, source observations, geometry and descriptors. Each call's support is its source observation plus distinct same-condition observations linked directly by saved strict anchor-overlap/same-direction matches at that grid. Multiple matched calls from one observation count once; cross-condition and transitive-only links do not contribute. The denominator includes every condition observation in metadata. Missing call tables are identified separately. Absent support products leave support counts `NA` and unverified; a validated empty match table gives source-only direct support. Anchor/body labels show counts, while blue indicates at least two supporting observations, dark gray source-only support and light gray unverified support. Native rectangles remain separate even when geometries overlap; no strengths are averaged or consensus coordinates created. The three original pairwise support products remain unchanged sidecars.

For local HiGlass conversion on Linux, create and activate the optional environment:

```bash
conda env create -f envs/browser.yml
conda activate microc2compare-browser
python scripts/export_browser.py --run /path/to/current/run --higlass
```

For conversion within the workflow, register this optional prefix with `BROWSER_ENV_PREFIX="/path/to/microc2compare-browser"` and enable `EXPORT_HIGLASS`. In locked environment mode, `BROWSER_LOCK` must name a tested platform lock with its recorded inventory. An explicitly configured prefix must pass the converter API and a small decoded bigWig smoke test during preflight. Without a browser prefix, the core environment is used and absent converters remain explicit inventory statuses. The table exporter itself does not install tools or activate environments.

The recipe pins UCSC `bedGraphToBigWig` 482, Clodius 0.22.2 and pyBigWig 0.3.26. Clodius 0.22.2 requires Python 3.12 or newer; pyBigWig 0.3.26 supports Python 3.9 or newer and provides Linux wheels. This is an optional conversion recipe, not a tested HiGlass server/client lock. Those requirements are recorded by the primary [Clodius package metadata](https://pypi.org/project/clodius/0.22.2/), [pyBigWig package metadata](https://pypi.org/project/pyBigWig/0.3.26/) and [Bioconda UCSC recipe](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/ucsc-bedgraphtobigwig/meta.yaml).

The browser recipe uses SciPy 1.17.1 with Python 3.12 and NumPy 2.2.6. Its [tagged release metadata](https://github.com/scipy/scipy/blob/v1.17.1/pyproject.toml) requires Python 3.11 or newer and NumPy 1.26.4 or newer below 2.7, covering these pins. The earlier SciPy 1.15.3 pin was unavailable to the actual Linux CI solve; successful solving and conversion still require the CI acceptance run.

The exporter uses the documented formats and command families:

```bash
bedGraphToBigWig signal.bedGraph profile.chromsizes.tsv signal.bw
clodius aggregate bedfile --chromsizes-filename profile.chromsizes.tsv --delimiter $'\t' --output-file intervals.beddb intervals.bed
clodius aggregate bedpe --chromsizes-filename profile.chromsizes.tsv --chr1-col 1 --from1-col 2 --to1-col 3 --chr2-col 4 --from2-col 5 --to2-col 6 --output-file rectangles.bed2ddb rectangles.bedpe
```

HiGlass uses `bigwig`/`vector`, `beddb`/`bedlike` and `bed2ddb`/`2d-rectangle-domains`. Clodius interval offsets depend on chromosome order. BigWig ingestion requires a matching chromosome-size coordinate system already registered on the same server; a remote chromosome track does not provide that registration. The bundle generates unique profile IDs from ordered names and lengths, converts interval databases separately for compatible profiles, and records structured commands in `higlass/ingestion_commands.json`. These requirements and commands follow the primary [HiGlass data-preparation documentation](https://docs.higlass.io/data_preparation.html).

The exporter never executes ingestion commands. Select the profile matching each displayed contact matrix before running a recipe against a server. Core vector, interval and rectangle tracks need no arc plugin; displaying loop arcs with the optional [HiGlass arcs plugin](https://github.com/higlass/higlass-arcs) requires loading that plugin and checking its compatibility with the actual client build.

Conversion states remain explicit. Missing executables produce `converter_missing` entries while preserving usable raw/UCSC products. A failing converter produces `export_failure.json` and no passed export audit. BigWig output gets a magic check and, when pyBigWig is installed, a decoded chromosome-header and interval/value comparison at the format's float32 precision. Clodius databases receive SQLite integrity checks; source-record decoding remains unverified. Receipts distinguish those checks from server ingestion and visual rendering. The installed HiGlass session, actual Linux conversion environment, plugin compatibility and rendering must be validated on the intended installation; no local table-export check proves them.

Validate the generic exporter and rolling sample scheduler with:

```bash
python -I -B -m pytest tests/test_browser_export.py tests/test_parallel_scheduler.py
```

The acceptance cases use a fresh unrelated two-versus-three cohort, hg38 and another reference identity, no prior exports, native loop geometry, disabled loop/annotation modules, multiple insulation windows, exact ADA IDs, optional stripes, corrupt inputs and explicit converter failures. Scheduler cases test free-slot refill, bounded workers, fail-fast cleanup, cached statuses and the older-Bash polling path. Actual descendant process-group cleanup has an additional Linux-only gate.

`tests/test_browser_pipeline.py` additionally runs the real CLI preflight, registration and browser stage plus the saved-results export action for five newly created synthetic matrices with two native chromosome orders. Scientific optional modules are disabled in this CLI integration; measurement export semantics are exercised by the acceptance fixtures. `tests/test_browser_runtime.py` checks the optional prefix contract with controlled converter doubles. The current Windows host has no pyBigWig, Clodius, UCSC converter or installed cooltools runtime; the initial local Linux recipe solve could not contact package channels. Hosted acceptance at `94f3e225af05` subsequently solved the recipe and passed 33 browser/scheduler tests, including actual bigWig, beddb and bed2ddb conversion and coordinate checks. Tested platform locks and intended-server ingestion/rendering remain separate installation checks.

The optional `.github/workflows/browser.yml` job installs the pinned Linux recipe and requires actual converters. `tests/test_browser_converters_actual.py` then exercises the exporter conversion path for bigWig, beddb and bed2ddb, decodes vector values, compares saved DB fields and exact cumulative chromosome offsets for two different native orders, and checks each Rtree entry contains the exact original geometry. The CI package cache is tied to the environment recipe; converted tracks are newly built under their coordinate profile IDs. The job retains package and conversion evidence. It does not ingest into a HiGlass server or test installed client rendering, and it has not been executed on this Windows host.

The CI job initializes its temporary directories in an early runner step and writes their paths to `GITHUB_ENV` for subsequent steps. GitHub permits the `runner` context in step settings, but not in job-level `env`; using `runner.temp` there prevents workflow validation. See the primary [context-availability table](https://docs.github.com/en/actions/reference/workflows-and-actions/contexts#context-availability) and [environment-file commands](https://docs.github.com/en/actions/reference/workflows-and-actions/workflow-commands#setting-an-environment-variable).
