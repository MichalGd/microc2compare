# Concurrent QC grids

QC processes each requested resolution independently. The grid set is the
distance-decay grid plus the enabled loop, Chromosight discovery and fixed-score
grids, with duplicates removed. A task may run several grids concurrently;
sample-level admission still reserves its full `THREADS_QC` and `MEMORY_MB_QC`.

Inner workers use separate spawned processes and open their own read-only
matrix handles. The worker count cannot exceed the number of grids,
`THREADS_QC`, `TOTAL_CPUS`, the effective CPU allocation, matrix/I/O limits, or
the task's RAM reservation. The parent reserves 512 MB. Each worker's estimated
reservation starts at 256 MB, adds 128 bytes per pixel in the actual streamed
chunk (`min(PIXEL_CHUNKSIZE, nnz)`) and 1,000 bytes per native matrix bin. Native
bin and pixel counts come from read-only Cooler headers, including native bins
outside the selected reference views. The largest requested grid estimate
determines admission. The estimate has no artificial upper cap; if the budget
cannot admit a worker, QC fails before calculation. `QC_WORKER_MEMORY_MB=0`
selects this estimate; a positive value explicitly overrides the reservation
for a reviewed workload. These are scheduling estimates, not operating-system
RSS limits; record actual peak memory and preserve headroom on a new dataset.
Native numerical-library thread counts are set to one before imports in spawned
workers.

The existing contact-summary JSON and expected-curve TSV paths, normalization,
masking, distance rules and numerical calculations are preserved. Expected
curves are produced from the fresh project's registered map, weights, common
mask and views. The run-owned cache checks those identities and the curve hash
before reuse. No previously completed cohort is required. Concurrent grids use
different cache directories, and failures do not publish a completed QC task.
HiCRep-only runs still emit their raw-input validation JSON without computing a
balanced curve.

The dedicated test compares actual two-grid serial and spawned calculations,
checks identical contact counts and expected tables, verifies unchanged input
checksums, and confirms that current cached curves are reused without rewriting
their files. This test establishes numerical equivalence on a fresh synthetic
custom assembly; it does not predict speedup or memory use on production maps.
