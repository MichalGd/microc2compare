# Validation and release gates

Local software verification and biological calibration answer different questions. The implementation must pass both the appropriate runtime tests and the biological/statistical gates before production inference. [Local verification](../validation/local-verification.md) records what actually ran; [environment solve evidence](../validation/environment-solves.json) records dependency resolution separately.

## Lightweight tests

Maintainers/developers may run these checks after selecting a testing prefix with the required dependencies (including pytest) and completing the dedicated-session setup in the installation guide. A cloned upstream prefix may not include pytest; do not install test-only packages into a published environment. Ordinary-user acceptance uses the smoke/full fixtures instead. Run from the selected source directory and keep test artifacts in your workspace:

```bash
cd "$m2c_src"
"$CONDA_EXE" run --no-capture-output -p "$CORE_ENV_PREFIX" python -m pytest -q \
  --basetemp "$m2c_work/validation/pytest-20261007" \
  -o cache_dir="$m2c_work/cache/pytest"
bash scripts/microc2compare.sh --help
find scripts -name '*.sh' -print0 | xargs -0 -n1 bash -n
```

Tests cover configuration/CSV quoting, resource admission, formula/rank/df rules, exact raw and inherited-weight queries, expected denominators including eligible zero pixels, zero versus NA, deterministic nontransitive consensus, missing weights/resolutions without repair, child failure and read-only integration. Additional restart/derived tests check output corruption, imported pool pixel sums, cleanup identity and seeded thinning. Scientific fixtures are generated before interception of balancing/construction APIs, so the guard exercises the workflow rather than preventing fixture setup.

`tests/test_handoff.py` uses original small Cooler files in the actual upstream layout to check merged/single-source selection, aliases and mixed-reference exports, run/sheet membership conflicts, chromosome mismatch, missing/failed/unknown weight evidence, integer storage, preserved identity gaps, CSV quoting, safe draft destinations and unchanged upstream hashes. A read guard proves preview does not read contact-pixel values. These tests validate input preparation, not biological independence or downstream calibration.

Default-path tests intercept balance/merge/zoomify/coarsen/create calls as errors after fixture creation. They compare input SHA256 values and check that no `.cool`/`.mcool` appears in outputs. The Linux scientific test additionally runs actual compartments/insulation/dots/aggregates against the virtual-weight adapter. R tests require explicitly available R packages and otherwise report a skip; a skip is not a passing R validation.

## Linux integration and CI

`.github/workflows/ci.yml` defines an Ubuntu core/HiCRep job and a separate R numerical adapter job. The core job checks Bash, Python/schema helpers, actual scientific adapters, the full descriptive demonstration, read-only input identities and resume. CI configuration existing locally does not mean GitHub executed it. Review uploaded logs and numerical results after publication. The isolated R job checks effect direction, exposure scaling, NA preservation and family-wide adjustment; it does not substitute for the ≥200-dataset calibration.

Interruption tests on Linux must verify that child groups exit, locks release, incomplete tasks have no success record and resume completes with unchanged imports. Inspect process-tree/RSS behavior on the actual server as part of the pilot, since OS/process tooling differs from Windows.

## Statistical calibration

[The prespecified protocol](../validation/calibration-protocol.md) and `scripts/calibrate.py` run the shipped R models on independent simulation datasets. Each supported design/module receives its own source/environment/recipe key. No passing artifact is supplied. Missing/mismatched evidence blocks formal/preliminary inference unless an explicit uncalibrated development run is requested and labeled. Add biological negative-control runs using genuine independent units; label shuffling/contact thinning cannot manufacture replication.

Inspect effect bias, interval coverage, dispersion/residual distributions, spatial dependence and family size alongside FDP/power. The v0.1 approximate loop uncertainty needs coverage review. Review weak effects, altered distance decay and broad redistribution against the declared relative-allocation truth. Record protocol revisions before rerunning; do not tune filters after seeing failure and reuse the old key.

## Public-data gate

`validation/datasets.tsv` identifies candidate papers/accessions without claiming downloads or completed analyses. Before using a candidate, register exact files/accessions, license/redistribution terms, raw and matrix hashes, assembly/reference assets, biological versus technical membership, batches/subjects, filters/balancing, independent feature provenance and expected biological comparisons. Reprocess upstream when deduplication/balance/axes are incompatible; downstream cannot repair them silently.

Compare adapter results with primary tool behavior at representative 10/100-kb resolutions, evaluate normalization/masking effects and biological negative controls, and review reproducibility against immutable source/environment identities. Candidate datasets are not automatically valid negative controls or independent feature sets.

## Release checklist

1. Pass Linux core/HiCRep/R numerical tests and full descriptive fixture; verify default no-balancing/no-copy behavior and resume/interruption.
2. Export/recreate verified environment locks and hashed pip artifacts; reject drift.
3. Pass prespecified simulation suites for every intended design/module; review uncertainty and biological negative controls.
4. Complete public-data/reference/membership checks and target-server resource pilot.
5. Review scientific outputs, document supported limits and decide project license/owner before publishing.

A local synthetic success establishes implementation behavior only. Deployment, publication and production readiness remain distinct decisions.
