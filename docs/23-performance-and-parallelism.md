# Performance and parallel execution

Version 0.3.3 optimizes repeated work in the 0.3.2 analysis path. These changes preserve chromosome scope, masks, balancing weights, candidate definitions, exclusion rules, event classifications and equal biological-replicate weighting. Imported matrices are opened for reading. Performance changes do not relax scientific eligibility criteria.

## What runs in parallel

| Stage | Actual execution | 0.3.3 improvement / remaining limitation |
|---|---|---|
| Preflight/register | Predominantly serial validation and I/O | Identity checks remain necessary; more CPU reservations do not accelerate file hashing |
| QC / HiCRep | QC sample jobs can overlap; streaming pixel calculations and HiCRep pair/view loops are serial within a task | Already completed in the historical pilot; chromosome/mask-specific results cannot be silently reused across different runs |
| Compartments / insulation | Bounded sample jobs and managed scientific-tool workers | The structural pilot completed these promptly; retain current measured settings |
| Mustache loops | Bounded sample jobs; `THREADS_LOOPS` drives a tile thread pool within each view | Remove an unused whole-map expected-contact calculation; native kernels and memory/I/O constrain scaling |
| Consensus | Serial deterministic grouping | Restrict same-resolution loop comparisons to active nearby groups; index usable boundary-derived intervals instead of repeatedly scanning the entire bin table |
| TADs | `THREADS_TADS` drives a **process pool across samples**; condition matching and export remain deterministic in the parent | Reuse sorted view/window arrays, binary searches and prefix counts across all windows and thresholds; index coverage and rejected-interval queries |
| Quantification | Bounded pixel streaming; vectorized matching, serial per sample/grid | Quantify both candidate sets in one scan; remove the Python loop over every stored pixel. `THREADS_QUANTIFY` remains an admission reservation |
| Chromosight | Sample jobs overlap; candidate windows are serial per sample/grid | Reuse expected-distance arrays; Cooler window reads and kernel calculations still cost time |
| Annotation | Shared anchors deduplicated; vectorized interval overlaps and indexed nearest TSS; final table formatting serial | Keep one full GTF/cCRE parse per stage; avoid repeated per-sample annotation or speculative extra workers |
| APA / ADA | Samples and feature windows serial within the aggregates task | Reuse expected-distance arrays. Cooler reads, ADA rescaling and repeated expected pixel scans across modules remain possible bottlenecks |
| Saddle | Serial sample/view accumulation of bounded matrices and partition groups | Uses the same reusable expected lookup; mandatory replicate and condition arrays and report plots are retained |
| Differential / report | R model fitting or descriptive summaries, then report rendering | Default single-worker R; biological calibration is separate from speed optimization. Plot rendering is serial |
| Pool / depth sensitivity | Explicit optional derived-map work | Disabled in the current fine profile; enabling these increases I/O and storage, rather than biological replication |

Execution CPU reservations are not utilization guarantees. OpenMP/BLAS/NumExpr workers remain fixed to one to avoid multiplying native threads inside workflow pools. The workflow stays within the configured allocation; it does not automatically occupy every server core.

## TAD settings and outputs

`THREADS_TADS` defaults to **2** and `MEMORY_MB_TADS` to **8000 decimal MB**. Existing project files inherit these values automatically. No new required user-supplied metadata or shell initialization is introduced. Optional overrides belong in the project configuration. Eight TAD workers are a reasonable acceptance trial for the reviewed 24-CPU, 160-GB allocation and eight observations:

```bash
THREADS_TADS="8"
MEMORY_MB_TADS="16000"
```

The actual worker count is bounded by requested workers, sample count, total CPUs, affinity/cgroup CPU allocation and estimated memory headroom. Parent tables, child copies and indexes are included in a conservative admission estimate. This is not an OS-enforced RSS ceiling. Each sample task computes all configured windows and thresholds; workers write no result files. The parent preserves sample/condition/event ordering.

`03_features/tads/performance.json` records actual workers, sample/window/threshold counts and domain-calling/total compute seconds. Existing scientific products and their meanings remain unchanged. `/usr/bin/time` and stage wall-time JSON supply complementary measurements; per-process maximum RSS does not establish aggregate process-tree peak RAM.

## Measured checks and practical tuning

`scripts/benchmark_performance.py --out NEW_JSON` compares the optimized algorithms with frozen 0.3.2 implementations on bounded synthetic data. It measures domain calling, event coverage, exact raw quantification and expected-window arithmetic, with numerical agreement checks. The quantification/window timings exclude HDF5 reads; they cannot predict full real-data stage duration. See `validation/performance-033-local.json` for the workstation measurements.

`scripts/validate_tad_optimization.py --baseline-run COMPLETED_PILOT --out NEW_DIRECTORY --workers 8` replays the accepted real-data insulation/mask inputs into a new private directory. It compares every prior TAD scientific product: coordinates, IDs, classes, reasons, counts and missing values must agree; floating scores allow only `1e-12` absolute/relative differences. It fingerprints inputs and prior outputs before and after. The completed pilot and source matrices are preserved.

`scripts/validate_performance_release.sh` combines that real-data TAD comparison with the existing actual Mustache/Chromosight/APA/saddle Linux demonstration. Passing this gate supports deploying the candidate, not interpreting biological differences or claiming calibrated statistical inference.

For the real fine-loop pilot, retain the accepted **24 CPUs / 160000 MB**, **two sample/matrix tasks**, and **four Mustache workers per sample** initially. Raising sample concurrency to four would reserve 16 CPUs and 128 GB with the current 32000-MB loop-task reservation. Try that only after the first measured sample tasks fit their reservation and storage throughput supports four readers. Increasing `THREADS_LOOPS` also multiplies the dense tile memory bound; do not remove the guard.

If APA/ADA or Chromosight then dominates the real pilot, measure candidate count, window size, task wall time, CPU, RSS and I/O before adding workers. Further candidates are bounded per-sample aggregate jobs, sample-local expected-table reuse with matching mask/weight identities, and locality-aware window reads. These are not implemented in 0.3.3. Arbitrary result reuse, unbounded window caches, and pooled biological samples would change risks or analysis semantics and are unnecessary for this optimization.

Versioned source updates invalidate completion signatures conservatively. Validate an old run with its original versioned command. Use a new run ID for the new scientific pilot; an ordinary `--resume` of an old-version run can recompute stages and is not a way to reuse accepted structural products across versions.
