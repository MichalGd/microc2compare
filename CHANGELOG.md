# Changelog

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
