# Running and restarting

Version 0.2 adds `--list-analyses` and `--analysis NAME[,NAME...]`, with automatic prerequisite switches and selected reports. See [the module catalogue](17-analysis-modules.md) for every module's inputs/parameters/outputs and named variants; [candidate upgrade](20-candidate-upgrade.md) describes installation beside the working release.

On a shared server, invoke the published release from any working directory with your own project configuration. Keep results/scratch outside the shared source and environments; pin the same release/profile/prefixes for resume. See [ordinary-user setup](16-shared-server-installation.md#3-start-as-an-ordinary-server-user). Only the maintainer installs shared dependencies. CPU/memory caps are per run, so concurrent accounts must follow the site's resource-allocation policy.

## Fixed stages

For a profile-bound shared installation, use `microc2compare -c /path/to/project.conf` directly, including `--analysis`, `--resume` and `--dry-run`. Its software profile loads automatically; [selection and precedence](25-automatic-software-profile.md) explain intentional overrides. Direct source-script examples below retain explicit software configuration behavior.

The runner follows this barrier-separated order:

```text
preflight → register → qc → pool → compartments → insulation → loops
          → consensus → tads → quantify → aggregates → differential
          → depth_sensitivity → report
```

Disabled optional stages record `skipped_config`. HiCRep runs as its own checked task at the QC barrier. Consensus/quantification run only when their enabled consumers require them. Core analysis does not depend on pool/depth tasks. Requested unsupported modules fail clearly. Small stage wrappers call `run_stage.sh`; `scientific.py`, `metadata.py` and `state.py` centralize their actual interfaces rather than duplicating parsing in each wrapper.

```bash
bash scripts/microc2compare.sh -c /path/to/project.conf --dry-run
bash scripts/preflight_check.sh -c /path/to/project.conf
bash scripts/microc2compare.sh -c /path/to/project.conf
```

Dry-run prints the selected enabled stages, threads, memory, concurrency and prerequisites without launching scientific computation. Preflight checks imported data/designs and environments read-only. Main execution additionally snapshots normalized settings and input identities. Initial SHA256 hashing of large imported matrices is necessary registration I/O; unchanged size/mtime reuses the cached SHA on ordinary resume.

## Completion and locking

One run root has an atomic `.run-lock` directory with an owner JSON containing PID, host and run ID. A second writer fails. A crash can leave a lock; check the owner process on the recorded host before explicitly removing that directory. There is no automatic stale-lock deletion.

Completion records bind stage/scope, relevant settings, source/schema/helper digests, imported/reference identities, environment inventory, current prerequisite signatures and required product sizes/SHA256. Product hashes are checked on resume. Registration arrays and TSV/JSON products are compact; stage commit points are written atomically only after products validate. A valid zero-call loop set has headers and a `success_empty` completion record. A failure never becomes a successful empty call set.

New pool HDF5 files publish from task-owned pending paths after raw sum/convergence validation. Imported HDF5 maps use read-only access. Failed intermediate files cannot satisfy a stage completion record and can be regenerated under the same explicit task request.

## Resume controls

```bash
bash scripts/microc2compare.sh -c project.conf --resume
bash scripts/microc2compare.sh -c project.conf --resume --until-stage loops
bash scripts/microc2compare.sh -c project.conf --resume --from-stage quantify
bash scripts/microc2compare.sh -c project.conf --resume --rerun-stage qc
```

`--resume` reuses only current valid products. Preflight reruns to refresh actual environment and reference checks. `--from-stage` requires valid prerequisites rather than manufacturing skipped records. `--rerun-stage` requires `--resume` and invalidates that stage plus its dependent commit points; it preserves inputs and existing products until replacements validate. Changing immutable assets/metadata should normally use a new run ID and clear provenance.

The source digest is deliberately conservative: some source/helper changes invalidate more work than strictly necessary. Skipping an expensive prerequisite after a source change requires review and a fresh validated record, never manually editing a `.done.json` signature.

## Failure and interruption

Every background child is collected with `wait`; a failed child fails its batch and the coordinator. Signals stop managed process groups and wait for children before lock cleanup. Logs are per stage/sample. A failed run writes `run_summary.json` with the failure state; successful earlier tasks can be reused. Linux process-group behavior still needs the Linux interruption test on the target runtime; Windows Git Bash is development verification only.

Task scratch is scoped to project/analysis/run/stage/sample, with `TMPDIR` pointing there. Success removes that owned directory when configured; failures retain it for diagnosis. Do not point scratch at imported-data directories. CPU/memory reservations govern admission, not a kernel-enforced per-task RSS ceiling; use a scheduler/cgroup for hard enforcement where required.

## Retention

By default there are no derived map tasks. When explicitly enabled, thinning maps are released per seed after P(s)/HiCRep/fixed-feature summaries; retention keeps only requested resolutions. New pooled maps are released after the report's consumers validate. Completion records mark discarded intermediates while retaining their scientific summaries. A stale pooled consumer regenerates only explicitly requested task-owned new pools. Imported maps and imported pools are never cleanup targets.

```bash
python3 scripts/validate_outputs.py --config-json /path/to/run/00_metadata/config.json
```

This verifies current enabled completion identities/hashes and HiCRep. For a deliberately partial run, supply `--until-stage NAME`. Use report status and per-comparison validation labels before interpretation; a completed software run does not establish calibrated biology.
