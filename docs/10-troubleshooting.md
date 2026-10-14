# Troubleshooting

For 0.3.0 missing fine grids, caller imports, dense-memory bounds, sparse windows, chromosome-scope migration and empty discoveries, follow [fine-resolution troubleshooting](28-fine-resolution-loops.md#resources-and-troubleshooting). Existing run scope cannot be changed in place; use its original versioned command or a new RUN_ID. Missing grids are never silently skipped, and masked outputs must not be interpreted as biological disappearance.

If acceptance stops at `compartments/C1.tsv differ` after both runs validated, measure oriented-score agreement before concluding that biology differs. Version 0.2.3 adds strict numerical/per-view raw-sign comparison; use [the correction guide](23-compartment-agreement.md). Meaningful score, coordinate, mask, component, eigenvalue and phasing disagreements still fail.

| Failure | Action |
|---|---|
| Software profile unreadable, failed, or pinned hash changed | Check the announced versioned launcher/profile pair. Keep published profiles immutable; select an intentional complete `--software-profile` override or have the maintainer publish a new pair. See [automatic profiles](25-automatic-software-profile.md) |
| Installation works only for its owner | Follow [shared access tests](16-shared-server-installation.md#4-accept-the-installation-from-a-second-account); check parent traversal, ACLs and all executable/prefix paths |
| Shared checkout/demo/config is not writable | Create the configuration and demo in your own workspace; source code/environments are deliberately maintainer-controlled |
| Several accounts exhaust CPU/RAM/I/O | Coordinate jobs with the site's scheduler/resource policy; per-run limits and run locks do not impose a server-wide budget |
| Upstream helper reports blocked/needs_review | Read its inventory and [handoff guide](13-upstream-handoff.md); complete actual identities/evidence, resolve missing or mismatched maps, then run normal preflight |
| Helper export refuses an existing directory or mixed references | Choose a new destination outside upstream results and export one reference with `--reference-id`; existing project files are never overwritten |
| Missing/invalid stored weight or explicitly failed convergence | Request corrected upstream output; downstream never repairs/rebalances it |
| Legacy convergence metadata absent | Obtain evidence upstream or explicitly document `finite_positive` legacy policy and its limitation |
| Missing requested resolution | Produce it upstream and register the corrected map; do not coarsen imported samples during core analysis |
| Reference/bin/filter mismatch | Correct assembly, chromosome order, view alignment, provenance or upstream processing before comparison |
| Split lanes separately deduplicated | Reprocess jointly per library upstream; matrix summation cannot remove cross-lane duplicates |
| Formula rank deficient or df <2 | Review batches/subjects, collect independent units, or explicitly request descriptive analysis; do not hide confounding |
| No independent frozen loop universe | Supply verifiable independent discovery or use descriptive same-data features |
| Missing/mismatched calibration | Execute/review the exact supported suite and register hashed passing evidence; development exception stays labeled |
| Core/HiCRep/R prefix check fails | Select a new owned compatible prefix and explicitly bootstrap; never update shared upstream prefixes during analysis |
| Bootstrap state checks raise `ModuleNotFoundError: No module named 'cooler'` through `state.py`/`matrix.py` | Version 0.2.0 accidentally imports a scientific module while computing registration signatures. Install [the 0.2.1 source correction](21-bootstrap-state-correction.md) and resume the existing fixture with unchanged environments |
| Locked mode missing inventory/drift | Recreate the reviewed lock and verify package inventory; no unpinned fallback |
| Configuration typo | Compare with the exhaustive variable reference; unknown uppercase assignments in the project file fail lint |
| Task exceeds configured budgets | Lower task workers/reservations only after measurement, or raise the verified allocation budget |
| Affinity/cgroup/memory allocation smaller than preset | Lower the global budget and every enabled task's reservation; logical CPUs do not override quotas |
| Dense view exceeds limit | Use suitable registered views, fewer workers or explicitly measured memory; keep the bound meaningful |
| Free disk floor fails | Point results/scratch at a suitable owned filesystem and budget optional derived data |
| Run root exists | Use `--resume` or a new stable run ID; avoid overwriting historical runs |
| `.run-lock` exists | Inspect owner PID/host liveness before explicit stale-lock removal |
| Completion stale or output corrupted | Resume prerequisites or use `--rerun-stage`; never edit signatures manually |
| Compartment export raises `TypeError: 'DataFrame' object is not callable` at `track(vectors, ...)` | The initial 7 October source has a phasing-variable name collision. Use the corrected 8 October source release and `--resume` with the existing project; no package change or imported-map rebuilding is required |
| Loop annotation raises `AttributeError: 'BinSelector' object has no attribute 'loc'` | Use the source release with the masked bin-selector correction, then resume the same project from preflight. This is an adapter type mismatch; do not remove the mask, change stored weights or reinstall packages to work around it |
| Comparison parameter writing raises `TypeError: Object of type int64 is not JSON serializable` | Use the source release with the design JSON correction, then resume the same project from preflight. The residual degrees-of-freedom value needs a standard integer; no R installation or package replacement is required for this descriptive fixture |
| Valid empty dots | Expect a header-only call table and `success_empty`; inspect coverage/depth rather than inventing calls |
| R loop validation raises `vapply ... result is length 0` | The old loop uncertainty adapter reads a missing edgeR field. Install the [0.2.4 source correction](24-r-loop-compatibility.md) and repeat R acceptance in a new private audit directory; do not change the RNA environment or skip uncertainty validation |
| Imported pool sum/member mismatch | Request exact upstream pool/provenance; no automatic repair |
| Optional new-map balance fails convergence | Inspect task-owned map coverage/recipe and rerun only its explicit task; imports remain untouched |

Start with `logs/STAGE.SCOPE.log`, `00_metadata/config.json` and the failing preflight message. The coordinator propagates every child's failure; a small exit-zero wrapper around a failed scientific command would invalidate the contract. Numerical warnings can indicate sparse/undefined measurements; inspect NA/status fields rather than suppressing them globally.

Paths with spaces must be quoted in Bash configurations/commands and quoted normally in CSV. The trusted `.conf` can execute Bash commands, so use only maintained project configurations. Linux does not resolve Windows drive letters; translate paths when transferring the project. Git Bash verification does not enable WSL or establish Linux R/cooltools support.

No generated result/plot is a substitute for input identity, biological independence or calibration evidence. Keep original upstream files immutable and preserve failed logs/scratch until diagnosis is complete.

When resuming with corrected source, its changed checksum can invalidate earlier completion records and rerun their tasks. Keep the same project configuration and generated inputs, start from preflight, and let the runner check prerequisites. Do not bypass checksum validation or force reuse of records made by different code. A new source release can reuse accepted software environments without reinstalling them.
