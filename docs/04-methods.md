# Scientific methods

For 0.3.0 Mustache/Chromosight adaptations, complete multi-resolution measurement, original-call preservation and projected/collapsed pixels, see [fine-loop methods](28-fine-resolution-loops.md). Independently configured insulation/TAD/ADA and optional compartment/saddle methods remain available. Autosome selection defaults on; stored upstream weights retain their original balancing scope. Reports separate the requested distance range from Mustache's effective near-diagonal discovery limit.

## Quantities and analysis domain

Raw counts are nonnegative integer unique-contact pixels. Balanced values multiply those counts by validated stored upstream weights; divisive weights are inverted in memory. Expected values are average contacts per eligible pixel at a view-specific genomic separation. Observed/expected (O/E), log2 insulation, boundary prominence and caller enrichment are distinct quantities. None substitutes for raw counts in the loop differential model.

The common domain combines registered nonoverlapping views, the blacklist and valid inherited weights across samples. Masked weights are sidecar arrays exposed through an original read-only Cooler API adapter. The input HDF5 files remain unchanged. Restricting the domain inherits upstream normalization; it does not balance the restricted submatrix. Missing values remain NA; an eligible observed zero remains zero.

The default core path does not create `.cool`/`.mcool` products. Dense regional operations are bounded by `MAX_DENSE_BINS`; compartment views additionally have a provisional concurrent-array memory estimate. No genome-wide balanced or O/E matrix is exported.

An explicit native raw-only HiCRep request sets `RUN_QC=false`, `RUN_HICREP=true` and disables all balanced feature/aggregate modules. The QC barrier then validates raw HiCRep inputs without producing balanced P(s); missing stored weights do not trigger repair or block that raw-only contract. Default `RUN_QC=true` includes balanced P(s) and therefore requires valid stored weights. The generator's `--raw-only` convenience flag disables cooltools features but leaves this default weighted QC enabled; it is different from a no-weight HiCRep request.

## Contact QC and P(s)

Counts are streamed in `PIXEL_CHUNKSIZE` blocks. Symmetric-upper storage is counted once; square storage uses one unordered upper representative and assumes true symmetry. Record this storage convention in upstream provenance; asymmetric square matrices are unsupported. Explicit new pooling requires symmetric-upper members; square-storage projects must supply validated upstream pools. Raw cis/trans totals differ from the scalar loop-model exposure, which includes only eligible same-view cis pixels in the inclusive 50-kb–2-Mb range.

Each P(s) row contains a view pair, distance, eligible-pixel denominator `n_valid`, raw/balanced sums and raw/balanced averages. An FFT of the eligible-bin indicator computes the denominator, including observable zero pixels. Distances below the ignored-diagonal setting have undefined averages. View averages must be combined by summing numerators/denominators, not by unweighted averaging. Current v0.1 calculates compact matching tables; it has no upstream expected-table import setting and never assumes a foreign expected table is compatible.

## HiCRep

The isolated HiCRep environment uses the package's sparse smoothing, stratum correlation and variance-stabilized weighting functions per common view. Input is raw counts, masked before smoothing; bins whose smoothing footprint touches excluded bins are excluded afterward. Counts are normalized within the requested view/distance range. Pair scores combine valid view SCCs using stratum weights. This avoids duplicate genome-wide contact files. It is an explicitly documented common-domain adapter and needs real-data comparison with the upstream reference implementation. Sparse/undefined scores remain NA, not zero. Pairwise concordance is descriptive and cannot establish biological replication.

## Compartments and saddles

At 100 kb, cooltools `eigs_cis` computes the first three components per view. The component with the greatest absolute Pearson correlation with the common immutable phasing track is selected, oriented positively, and accepted only when absolute correlation is at least `MIN_PHASING_CORRELATION`. Eigenvalues, component choice, correlation and ambiguous-view status are exported. An ambiguous view has NA oriented scores; there is no fabricated orientation.

The default descriptive saddle partition uses each sample's oriented percentile ranks, average ranks for ties, and an equal-weight, condition-blind mean within each view. Only bins finite in every sample participate. Fifty fixed groups use deterministic bin-ID tie breaking; views with fewer than fifty usable bins are flagged. Partition membership and its hash are frozen before per-replicate saddle aggregation. This proposed partition does not claim equivalence to eigendecomposition of a pooled map. Alternatives explicitly require an external exact-grid reference track or a validated all-sample pool.

Saddles aggregate view-matched O/E and retain eligible-pixel counts. AA/BB/AB summaries and compartment changes are descriptive. Domain-sized dense slices respect the configured bound. No dcHiC or other formal compartment model is substituted.

## Insulation, boundaries and domains

The 0.2 [TAD module](18-tad-analysis.md) adds per-replicate and condition-supported insulation-defined domains and descriptive event groups. Insulation retains cooltools' Li threshold by default; Otsu/fixed thresholds and fixed-threshold TAD sensitivity are explicit settings. Domain/event coverage, matching, replicate disagreement and length changes are recorded. A geometry event has no inherited bin/caller FDR.

At the fixed 10-kb grid, cooltools insulation calculates 50-, 100- and 250-kb windows by default. The 100-kb window is the primary inferential family; additional windows are descriptive sensitivity products. Common masked weights, matched views and ignored diagonals are passed to the API. Scores/boundaries also require the prespecified valid-pixel fraction relative to the theoretical full diamond; edge/mask deficits remain missing.

Boundary prominence is present only where supplied by the algorithm. Non-minima are not assigned zero strength. Deterministic boundary consensus clusters locations within one bin, per view/window, without transitive drift; domains are adjacent consensus boundaries with no excluded-bin bridge. Domain strength is mean intra-domain O/E with an explicit eligible-pixel denominator. Domains and boundary-location changes remain descriptive; a merged interval does not acquire regional FDR from individual bins.

## Dots, consensus and raw counts

cooltools `dots` runs at 10 kb with view-matched expected tables and the masked-weight adapter. Retained cis anchors span 50 kb–2 Mb. Caller enrichment significance describes discovery; it is distinct from biological differential P values. Optional 5-kb calling requires the explicit reviewed resolution gate and existing valid 5-kb input groups. Finer calling is rejected/deferred.

Anchors are canonical single-bin BEDPE coordinates, ordered deterministically and deduplicated. Consensus uses complete-linkage tolerance of one bin on both anchors; chained neighboring calls cannot expand a cluster beyond that tolerance. The medoid, coordinate-derived feature ID and biological/source support are retained. Same-data consensus is descriptive even when discovered from a condition-blind pool.

Quantification streams original raw pixels and measures one anchor-pair pixel per feature. Eligibility requires the common mask, same view/chromosome and configured distance range. Missing/excluded features carry reasons and NA; eligible absent stored pixels are zero. The exposure sums each unique eligible unordered raw cis contact across all views, not the subset of feature counts.

## APA and ADA

APA averages matched O/E windows with a fixed 100-kb flank around consensus loop anchors, preserving missing pixels and contribution counts. ADA rescales boundary-derived domains to fifty bins and averages eligible O/E; domain strength separately records its pixel denominator. Per-replicate arrays and summary TSVs retain the descriptive status. Empty feature families yield schema-valid products with NA aggregate measurements.

## Optional derived products

Pool import streams an exact sum validation and reads existing pool weights. Explicit creation makes only the required analysis resolutions, checks sum conservation, and balances only genuinely changed, task-owned counts after its separate opt-in. A single-member pool aliases its sample. Imported files and weights never change.

Depth sensitivity is off by default. When enabled, seeded NumPy PCG64 binomial thinning samples the finest requested eligible cis grid, then derives only requested integer coarsenings from that same draw. Target depth is the minimum native eligible sample depth; realized totals are recorded because binomial thinning matches expected, not exact, depth. Every seed gets raw P(s) and, when enabled, HiCRep. `raw_only` creates no new weights. `balanced_features` requires explicit new-map balancing, matched expected values and fixed-loop/primary-insulation summaries. Identical count draws reuse inherited weights without balancing a copy. No derived pool/thinned map enters a replicate-level inferential model.

Task-owned derived maps are released after their consumers succeed, unless retention is requested. New pool records survive cleanup with explicit discarded-intermediate annotations; a stale consumer regenerates those maps only under the original explicit creation/balancing settings. Imported pools are never removed.

Primary API sources: [Cooler](https://cooler.readthedocs.io/en/latest/api.html), [cooltools](https://cooltools.readthedocs.io/en/latest/), [HiCRep](https://github.com/dejunlin/hicrep). Exact numerical adapters and synthetic tests are in `src/microc2compare/` and `tests/`; production equivalence still requires the separate validation gate.
