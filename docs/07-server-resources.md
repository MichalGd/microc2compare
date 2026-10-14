# Server inventory and resource tuning

The reported server has 500 GB RAM, four Xeon processors, 70 physical cores and 140 logical threads. These values were not verified locally. The Windows workstation is a different machine. Use the allocated resources and measured process behavior; logical threads are not additional physical cores.

For the [shared installation](16-shared-server-installation.md), budgets and preflight checks apply per run, not across all accounts. Two users cannot each claim the full machine. Use the existing scheduler with matching CPU/memory requests; without one, coordinate a central submission policy and run only one heavy pilot at a time until it is established. Run locks protect one output tree and do not reserve global resources. Keep each account's results/scratch/cache on its own approved writable paths.

```bash
m2c_src="/opt/conda_envs/microc2compare/releases/shared-design-json-fix-20261008/microc2compare-0.1.0"
m2c_project="$HOME/microc2compare-work/projects/pilot"
umask 077
mkdir -p "$m2c_project"
bash "$m2c_src/scripts/server_inventory.sh" > "$m2c_project/server-inventory.txt"
```

The script records OS, CPU/NUMA topology, memory, filesystems, limits, affinity and cgroup information without changing the server. Review scheduler allocations, storage mount type, free disk, inode counts, tmp space and shared workload. Missing optional inventory commands are reported rather than installed automatically. No NUMA tuning, service change or kernel setting is applied.

## Presets

| Preset | CPU budget | Memory MB | Jobs | Samples/matrix/I/O |
|---|---:|---:|---:|---|
| demo | 4 | 16000 | 4 | 1/1/1 |
| pilot | 24 | 160000 | 8 | 2/2/1 |
| production candidate | 60 | 380000 | 12 | 4/4/2 |

Memory is decimal MB. The production budget leaves headroom relative to the reported capacity but remains provisional. Source the appropriate preset explicitly in the trusted project configuration and override individual `THREADS_*`/`MEMORY_MB_*` values afterward. Demo tasks reserve one thread/8000 MB, including optional balancing, so they fit the demo budget. A smaller CI/laptop allocation needs smaller task reservations as well as a smaller global budget.

Batch concurrency is the minimum of total-job cap, sample cap, `floor(TOTAL_CPUS/task_threads)`, `floor(MAX_MEMORY_MB/task_memory)`, and applicable matrix/I/O caps. Stages use homogeneous bounded batches with barriers. New-map tasks conservatively reserve the larger of their own and balancing requirements. A task larger than the global budget fails before launch. Reservations are admission policies; they cannot force a tool to honor an RSS ceiling.

Preflight rejects CPU budgets exceeding process affinity/cgroup-v2 quotas and memory budgets exceeding current physical/cgroup available memory. It checks every visible cgroup ancestor. Cgroup v1 is explicitly unsupported; use a reviewed v2 allocation or unrestricted host. It also requires the configured minimum free disk near results/scratch. That floor does not estimate all optional derived-map space: calculate map/resolution/seed size from pilot measurements before enabling those tasks.

Nested OpenMP/BLAS/NumExpr threads are fixed to one. cooltools dots/insulation receive their task worker count; eigendecomposition and new-map balancing use managed process pools. TADs use a bounded spawned process pool across samples, with THREADS_TADS and MEMORY_MB_TADS controlling admission. Quantification uses vectorized bounded pixel chunks and measures candidate sets together. QC pixel streaming, quantification sample loops, aggregates and HiCRep still contain serial work, so their configured thread values are not promises of full CPU utilization. See the [stage-specific performance review](23-performance-and-parallelism.md). Thinning is serial by seed/sample to bound intermediates. R models use the reserved single worker by default.

## Pilot sequence

First validate synthetic Linux behavior, then run one representative biological analysis set at 10/100 kb with pools/thinning disabled. Confirm upstream assembly, masks, weights and reference identities. Use the pilot preset only when the actual allocation permits it. Start with two sample tasks and one I/O-heavy writer. Increase workers and matrix concurrency only after measuring wall time, CPU use, RSS, I/O, free disk and feature sizes.

```bash
bash "$m2c_src/scripts/bootstrap_environments.sh" -c "$m2c_project/config/project.conf" --check-only
bash "$m2c_src/scripts/preflight_check.sh" -c "$m2c_project/config/project.conf"
bash "$m2c_src/scripts/microc2compare.sh" -c "$m2c_project/config/project.conf" --dry-run
bash "$m2c_src/scripts/microc2compare.sh" -c "$m2c_project/config/project.conf" --until-stage qc
# Inspect QC and benchmark evidence, then resume the approved descriptive stages.
bash "$m2c_src/scripts/microc2compare.sh" -c "$m2c_project/config/project.conf" --resume
```

Each task writes wall-time JSON. On Linux, `/usr/bin/time -v` additionally records maximum RSS, CPU and filesystem I/O in `benchmarks/*.time.txt` when available. Its maximum RSS is not necessarily aggregate simultaneous process-tree RSS; measure the scheduler/cgroup peak or monitor all child processes for parallel tools. Windows wall-time results do not establish Linux resource requirements. No target-server benchmark is supplied.

For compartments, large views can require substantial dense memory. `MAX_DENSE_BINS` and the twelve-array concurrent estimate gate requests. Split biologically appropriate chromosome-arm views, reduce workers, or increase explicitly measured reservations; do not blindly remove the bound. APA/domain/HiCRep feature distributions and exact imported-pool validation can dominate runtime. Profile those before production scale. Optional pooling/thinning multiplies I/O and scratch requirements and does not increase biological sample size.
