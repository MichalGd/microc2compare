# Changelog

## 0.4.0 — general-purpose CPU, browser and stripe candidate, 13 October 2026

- Support new projects with their own assemblies, observation counts, conditions and registered matrix grids. Add generic `hsieh_8grid` and `short_range_4kb` presets, input capability reports and explicit exact-position retained-pairs preparation of genuinely new maps. Fine-grid admission remains project-specific; coarser pixels and weights are never promoted to authentic finer grids.
- Add independent CPU Chromosight discovery alongside fixed-candidate scoring. Retain caller identity, native coordinates, provenance and cross-resolution correspondence; multiple callers from one observation do not increase biological support. Document the current Mustache six-bin ideal occupancy floor and the limits of targeting 4 kb genomic separation. These are adapted profiles, not an exact historical Hsieh reproduction or a detection guarantee.
- Add bounded CPU/RAM admission, a rolling sample scheduler, per-grid QC workers, parallel quantification/aggregation, atomic unit checkpoints and validated expected-contact reuse. Vectorize expected-table construction while preserving the existing streaming sums and denominator definitions. Resource estimates require representative runtime measurements; they are not hard RSS limits.
- Reuse dependency validation and source hashes within one audit, with unchanged signatures and fresh checks on subsequent calls; prevent report validation from repeatedly traversing identical prerequisites.
- Add optional StripePy 1.3.0 discovery in a separate environment. Temporary raw cis maps respect common usable-bin segments; original matrices remain unchanged. Preserve native stripe anchors/bodies, descriptors, direction, coverage and descriptive cross-observation support without claiming biological significance.
- Add dated browser bundles for saved individual and condition results, including native loops, all saved insulation windows, boundaries/TADs, exact ADA domains/scalars, saved compartment components and stripes. Retain geometry, missingness, reference profiles and source checksums. Optional UCSC/HiGlass conversion has explicit status; export does not ingest into a server or certify rendering.
- Add explicit donor/culture metadata for descriptive observations and reject unsupported nested-culture inference. Preserve historical run identities; install the candidate alongside the accepted release and use new scientific run IDs.

Local numerical/regression and synthetic integration checks include unrelated projects, unequal condition sizes, multiple references, input checksum preservation, serial/parallel agreement, checkpoint recovery, Mustache geometry/band access and actual Chromosight APIs and isolated one/two-process StripePy execution. Public-checkout checks cover 372 passing cases across the full suite and its corrected-fixture rerun; 10 checks were skipped. The full run recorded 371 passes, 10 skips and one synthetic peak-localization fixture failure; the corrected fixture plus CPU caller/pairs tests passed 27 checks. These overlapping counts are not summed. Product sources were unchanged during the correction. Linux CI execution, optional converter/runtime acceptance, shared-server resource measurements, real-data short-range/stripe sensitivity, installed-browser rendering and formal biological calibration remain separate pending checks. Earlier release acceptance below does not establish 0.4.0 Linux acceptance. See [the upgrade guide](docs/35-general-purpose-upgrade.md), [CPU loops/pairs](docs/cpu_loop_profiles.md), [browser export](docs/browser_export.md) and [stripe detection](docs/34-stripe-detection.md).

## 0.3.3 — performance candidate, 10 October 2026

- Index TAD domain/coverage/disagreement queries and reuse arrays across thresholds; add bounded spawned sample workers with independent CPU/memory settings and execution metrics.
- Measure both loop candidate sets in one vectorized pixel scan while preserving integer counts, eligible exposures and exclusions.
- Restrict loop consensus to active nearby clusters and index boundary-derived domain masks.
- Avoid unused expected-contact computation for Mustache; reuse expected-distance arrays for Chromosight, APA, ADA and saddle windows.
- Add frozen 0.3.2 numerical comparisons, synthetic algorithm benchmarks and read-only replay acceptance against a completed structural pilot. Existing releases/matrices/configs remain intact; server acceptance is required before activation.

## 0.3.2 — project convenience, 10 October 2026

- Add shared-command synthetic QC demonstration and literal parameter variants.
- Add compact project initialization with installed scientific presets and explicit per-project overrides.
- Route upstream discovery/combination, reference/annotation registration, resolution inventory/preparation and completion validation through the same profile-bound `microc2compare -c project.conf` interface.
- Derive fixed setup and repeated helper arguments internally. Preserve draft evidence gaps, existing hash checks, biological units, imported map hashes and new-directory/no-overwrite behavior.
- Reuse registered prepared grids; generate companion configurations automatically. Scientific methods, required saddles and statistical gates remain unchanged.
- Replace obsolete ordinary-user clean-shell/profile/cache instructions with the shared launcher interface; retain installation/acceptance steps for maintainers.

## 0.3.1 — required saddle plots and final reports, 2026-10-10

Full practical/Hsieh-inspired profiles now require compartments, saddle aggregates and the final HTML report. Explicit `--analysis` selections remain independent. The fine-loop example enables 100-kb compartments and APA/ADA/saddle by default. Preflight validates the fixed phasing input before contact scans; a saddle with no eligible common oriented partition fails clearly. Full reports require current arrays for every biological replicate and condition, embed their PNG plots, link SVG/NPZ outputs, and show per-replicate segregation summaries. Condition plots retain equal biological-replicate weighting.

Extended Linux acceptance includes actual compartment analysis, nonempty saddle arrays, embedded replicate/condition plots, saddle-only execution and full/selected numerical agreement. The small synthetic fixture uses an analytically defined phasing signal and a 10-kb compartment grid to provide at least 50 bins; real-data defaults remain 100 kb. Existing core/HiCRep/R and optional caller dependencies are unchanged. Install this version alongside 0.3.0; target-server acceptance and real-data phasing verification remain pending.

## 0.3.0 — fine-resolution loops and autosome policy, 2026-10-10

Added pinned Mustache 1.3.3 discovery, complete 2/4-kb candidate quantification and adapted Chromosight 1.6.3 fixed-window scores. Cross-resolution matching preserves source calls, exposes ambiguous correspondence and prevents duplicate coarse-pixel tests. Added both-anchor GTF/cCRE annotation, per-grid APA/comparisons/volcanoes, equal-replicate condition aggregates, and APA/ADA/saddle PNG/SVG with numerical coverage and shared plot scales. ADA retains its own grid. Added read-only resolution inventory, explicit preparation of only missing integer coarsenings with new-map balancing, reviewed 1-kb/sub-kb admission and primary-source method metadata.

Autosomes are the default for mouse and human; X/Y are opt-in. Frozen chromosome/view identity prevents scope migration on resume. Launchers handle private permissions and scratch-local caches. Shared core/HiCRep/R are reused; optional caller software has a separate checked prefix. Historical runs require their original versioned command or a new RUN_ID. Scientific signatures/calibration recipes change. Source is a validation candidate; target Linux Mustache acceptance, server deployment, real-data pilot and statistical calibration remain distinct pending steps.

## 0.2.7 — legacy input manifests, 2026-10-09

The input helper now reads older sample/merge manifests without reference columns and the known headerless three-column technical membership format. Missing reference evidence remains explicit; known conflicting references and malformed records still fail. Added opt-in `--individual-matrices` for a sheet with exactly one source per declared biological unit: select each existing individual map while preserving historical pooled grouping/replicate conflicts and missing sample provenance for review. The default merged-input rules, subject/library checks, draft evidence gaps and no-overwrite behavior remain strict. Exported inventories and provenance sidecars retain old records and review notes. Scientific stages, their source fingerprints, software profiles and environments are unchanged from accepted 0.2.6. No real matrices or installed source are modified; target-server helper acceptance and real-data provenance review remain required.

## 0.2.6 — common chromosomes and X/Y, 2026-10-09

Added exact coordinate alignment for imported maps with different native chromosome orders or extra unselected scaffolds. Registered views select chromosomes, including X/Y when requested; absent selected chromosomes, unequal lengths/bin grids and incompatible count/filter/weight recipes still fail. Each native map retains its stored weights and raw counts. Feature tables use the common reference order; cooltools receives native-order views. Added scope/coverage audits and report explanations, explicit QC total/exposure scopes, common-grid optional pool creation and thinning, and missing assembly-metadata review in the input helper. Imported-pool verification still requires identical complete native member grids. Added numerical regressions and synthetic Linux scientific acceptance instructions. Existing environments, inputs, shared releases and active launcher are not modified by archive installation. Source signatures/calibration keys change; server acceptance and statistical calibration remain separate requirements.

## 0.2.5 — automatic shared software profiles, 2026-10-08

Added profile-bound versioned launchers, pinned profile hashes, automatic software-only loading after project configuration, explicit `--software-profile` overrides and recorded profile provenance. Ordinary users no longer need a fixed manual software `source` line. Project inputs, outputs, scientific parameters and resource budgets retain their own configuration values. Added actual launcher/configuration/QC regressions and user/maintainer upgrade instructions. Existing software environments, profiles, releases and scientific/R algorithms are unchanged; the common-helper change conservatively invalidates stage source fingerprints, so resume may recompute stages. Target-server interface acceptance remains required.

## 0.2.4 — edgeR loop uncertainty compatibility, 2026-10-08

Fixed the R loop adapter's missing `var.post` access with edgeR 4.4 by selecting `s2.post` when available. Validated scalar/per-feature dispersions, used the fitted working NB dispersion for approximate uncertainty and NB Pearson residuals, and added explicit diagnostic columns. The edgeR fit/test, filtering, exposure, effect and P/FDR calculations remain unchanged. Added actual-R comparisons with direct edgeR and analytical variance plus invalid-field/scaling regressions. Other descriptive scientific fingerprints remain compatible with 0.2.3; differential fingerprints and calibration keys change. Server acceptance and intended-design calibration remain required. Existing archives, shared environments and the active 0.2.3 command remain unchanged.

## 0.2.3 — numerical compartment acceptance, 2026-10-08

Replaced byte comparison of full/selected compartment output with a dependency-free numerical bundle check. It permits per-view raw eigenvector signs and small rounding differences, while preserving oriented score/browser direction, coordinates, missing patterns, component identity, eigenvalues and phasing consistency. The observed server oriented scores already agree within 1.6e-15. Scientific algorithms, models, environment/configuration contracts and scientific-signature source files remain unchanged, allowing compatible completed results to be reused. Added regression/CLI/Bash coordination tests and safe resume instructions; acceptance remains pending on the target server.

## 0.2.2 — differential loop volcano reports, 2026-10-08

Added report-only volcano plots from current completed biological loop models: edgeR log2 relative-allocation effects versus differential P values, coloured by the existing full-family BH FDR. PNG/SVG/point TSV exports and status/filtering counts preserve original results. Descriptive comparisons explain why no volcano is appropriate; preliminary/development labels remain explicit. Optional display thresholds/labels do not change statistical tests. Plot exports are linked from the HTML and documented with standalone-module/report refresh commands. This amended candidate includes the 0.2.1 bootstrap fix and existing-fixture resume support; prior archives/installations remain unchanged.

## 0.2.1 — bootstrap state correction, 2026-10-08

Moved the shared existing-resolution selection policy into the dependency-free analysis planner. Bootstrap prerequisite, completion, summary and output checks no longer import Cooler/NumPy/pandas through matrix.py. Scientific registration imports the same policy, preserving resolution/weight semantics. Regression tests exercise the real state/validator CLIs while scientific imports are explicitly blocked.

Added `--resume` to the acceptance helper: it validates original matrix checksums and the existing software-profile reference, retains configurations and resumes existing/fresh selected runs without regenerating input maps. Missing identity files, mismatched profiles and changed input files fail before configuration writes. The source correction needs no environment changes; install it as a separate release and resume the failed candidate fixture.

## 0.2.0 — modules and TAD validation candidate, 2026-10-08

Added named analysis selection, selective APA/ADA/saddle consumers and configuration variants with distinct run IDs. Module selection enables prerequisites and avoids unrelated matrix resolutions, phasing, comparisons and R requirements. Registration no longer inherits unrelated preflight parameter changes; compatible results still require current inputs, software and prerequisite signatures.

Added insulation-defined per-replicate/condition TADs, support/disagreement/exclusion tables, descriptive merge/split/growing/shrinking/appearing/disappearing/shift/strength/uncertainty events, fixed-threshold sensitivity and browser BEDs. Preserved Li boundary threshold as the default. Reports show selected settings/input identities, TAD summaries and bounded boundary overlays. Added the complete module/output catalogue, R acceptance guide and candidate upgrade steps with explicit Windows/Linux commands.

Added read-only R environment discovery/numerical acceptance, reviewed-existing calibration registration, verified non-overwriting release installation and explicit versioned launcher publication. The working server release/environments are not modified by local development. GitHub publication is deferred. New Linux scientific acceptance, ordinary-account execution, actual R/calibration and biological pilot evidence remain required.

## Unreleased — design JSON correction, 2026-10-08

Converted residual degrees of freedom to a standard Python integer when constructing the design summary so strict JSON parameter writing works in descriptive and inferential modes. Statistical formulas and JSON validation remain unchanged. Added descriptive loop/insulation output checks and design-serialization coverage for all inference modes. Updated installation examples to select the corrected source candidate while retaining the existing environment prefixes.

## Unreleased — masked bin-selector compatibility correction, 2026-10-08

Made the read-only masked bin adapter a Cooler `RangeSelector1D`, allowing the real `cooler.annotate()` call in cooltools dots to recognize and slice it correctly. Virtual analysis weights, including excluded-bin NaNs, remain in memory; no matrix copying, construction or balancing is added. Added annotation, bounded-query and selector regression checks, and removed an unnecessary HiCRep dependency from the core-only cooltools integration test.

## Unreleased — compartment export correction, 2026-10-08

Fixed a local phasing DataFrame hiding the browser-track writer in the compartment stage. The eigensolver, stored-weight handling and orientation rules are unchanged. Added a regression check for oriented, missing and ambiguous bins and all four exported products. Server inventory now avoids environment-enumeration commands that can clean a user's existing Conda registry.

## Unreleased — upstream handoff and documentation, 2026-10-07

Made shared all-user server installation the documented deployment model: maintainer-controlled source/environments and versioned software profile/command, separate user configurations/results/scratch/cache, non-maintainer acceptance and cross-user resource coordination. Moved public demo/quick-start writes outside the checkout. Bash helpers disable Python bytecode writes in shared software. Added `config/site.example.conf` and the shared-installation guide; actual server access remains to be tested.

Added a numbered server-installation walkthrough with explicit workstation/server commands, existing Conda/Python/prefix discovery, a software-only bootstrap configuration, selective isolated core/HiCRep fallback creation, deferred differential R, conservative resource budgets and smoke/full-descriptive fixtures before upstream piloting.

Added a read-only microc2tracks input-preparation helper with manifest-aware biological-unit matrix selection, resolution/count-storage/weight inspection, related-output inventory, explicit evidence gaps and non-overwriting draft exports. Added real Cooler fixture tests, a worked upstream sheet, README Mermaid diagram, beginner walkthrough/glossary and documentation index. Existing runner/configuration/output conventions and statistical gates remain unchanged. Linux/R/calibration/public-data/server gates are still pending.

## 0.1.0 — implementation candidate, 2026-10-07

Implemented a fixed shell-based stage runner, immutable upstream count/weight registration, resource admission, restart signatures, Python feature/QC/aggregate adapters, guarded optional pools/thinning, edgeR/limma adapters, HTML reports, schemas, examples and installation/interpretation documentation.

Local Windows read-only synthetic integration and lightweight tests were executed. Linux environments resolved; full Linux/R execution, exact tested locks, calibration, public-data validation and the server pilot remain explicit release gates. No server deployment, Git initialization or GitHub publication occurred.
