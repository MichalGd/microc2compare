# First real-data analysis

Use your ordinary Linux Bash login. After the maintainer installs and accepts the current profile-bound command, no Conda activation, clean-shell setup, manual profile source, software paths or runtime exports are required.

```bash
microc2compare --init "$HOME/my_microc_project/config/project.conf"
```

Edit the created configuration's project/genome, input locations, output/scratch directories and CPU/RAM budget. The fine profile defaults to 2/4-kb loops, 10-kb insulation/TAD/ADA, 100-kb compartments and APA/ADA/saddles in the final report. Use `--preset existing` during initialization for the older 10-kb profile. Existing finalized projects skip initialization and reuse their completed inputs.

Follow [simple project commands](32-simple-project-commands.md) to import upstream results if needed, complete genuine sample/provenance facts, reuse/register references and optional annotations, then audit requested resolutions. The helper inventories drafts; it does not establish biological independence or approve formal inference. [The upstream handoff guide](13-upstream-handoff.md) explains its evidence fields. Register reference hashes automatically rather than typing them; a common usable phasing track is required for full compartment/saddle analysis.

```bash
microc2compare -c /your/project/config/project.conf --audit-resolutions
microc2compare -c /your/project/config/project.conf --analysis qc,hicrep --dry-run
microc2compare -c /your/project/config/project.conf --analysis qc,hicrep
microc2compare -c /your/project/config/project.conf --validate
```

If a preparation operation prints a new sibling configuration, use that path. Inspect contact totals/retention, raw and balanced distance-decay curves, usable-bin coverage and biological-replicate concordance before the full run. All results are under `OUTDIR/PROJECT/ANALYSIS_SET_ID/RUN_ID`; open its `08_reports/index.html`. Undefined measurements stay missing, not zero. Choose a new run ID for a separate history, or use reviewed `--resume` behavior to extend a selected run. The [module guide](17-analysis-modules.md) explains separate analyses and prerequisites.

```bash
microc2compare -c /your/project/config/project.conf --resume
```

The full fine-loop run includes obligatory saddle plots in its HTML report. It needs an accepted caller environment, valid resolution-specific inputs and the documented scientific gates. Only an explicit `--prepare-resolutions` operation creates missing exact coarsenings and weights for those new maps; analysis never rebalances imported maps. Running through `nohup` can keep long jobs alive after SSH disconnection; see the project guide.

The synthetic demonstrations test software behavior, not biological calibration. The maintainer runs full/selected/resume acceptance from [the current upgrade guide](33-simple-project-upgrade.md) before activating a candidate. That one-time process uses installed paths and private demo destinations, rather than requiring every researcher to repeat installation setup.

## Read effects before significance

A loop can rise from 100 to 180 raw contacts while the sample's eligible contact total doubles. Its relative allocation then falls slightly: `180/(2L)` is smaller than `100/L`. The model measures relative allocation, not absolute interactions. Balanced values help feature analysis and visualization; raw counts supply the loop statistical model.

Loop-caller significance concerns contact enrichment in one map. It is not evidence that treatment changed the loop across biological replicates. Exploratory results have biological differential P/FDR missing. Formal results need the declared design, independent loop candidates (or fixed insulation grid), calibration and all other gates. See [statistics](05-statistics.md) and [worked output interpretation](09-outputs.md).

Browser files currently include bedGraph signals for compartments/insulation and interval/loop tables. Load them using the registered assembly and exact coordinates; use upstream `.hic` files for contact-map viewing. See the upstream browser guide for `.hic` publication. Automated downstream BigWig or track-hub publication is not claimed.

## Short glossary

| Term | Meaning here |
|---|---|
| Contact matrix | A table of contact counts between pairs of genomic intervals, stored sparsely |
| Cooler / `.cool` / `.mcool` | Contact-map format; `.mcool` stores several genomic bin sizes in one file |
| Resolution / bin | Width of the genomic intervals; 10 kb means 10,000 bases per interval |
| Balancing weight | Upstream bin correction factor; it does not replace a raw count |
| cis / trans | Contacts within one chromosome / between different chromosomes |
| P(s) | Average contacts as genomic separation increases; a contact-decay summary |
| HiCRep / SCC | A contact-map similarity method and its stratum-adjusted correlation score |
| Loop / dot | A focal enriched contact between two genomic anchors |
| BEDPE | A tabular representation of two genomic intervals, such as loop anchors |
| Insulation / boundary | A local contact-crossing score / a candidate separation between contact domains |
| Compartment | A broad genomic interaction pattern, oriented consistently using a reference signal |
| Consensus feature | A deterministic feature set assembled from discovery calls; its provenance controls interpretation |
| Pool | A combined discovery/visualization map; never an additional biological replicate |
| Exposure offset | The statistical model's declared contact-total denominator for each biological sample |
| FDR | False discovery rate: the target expected fraction of false discoveries among reported discoveries under the method's assumptions |
| APA / ADA / saddle | Aggregate contact summaries around loop anchors / domain features / compartment groups |
| Provenance | Records of the actual inputs, processing, versions, reference and parameters |

## If something fails

Start with the exact preflight message or `logs/STAGE.SCOPE.log`, then [troubleshooting](10-troubleshooting.md). Missing input weights/resolutions require upstream correction; missing biological information requires actual experimental records. Do not turn a failed task into a successful empty result or change evidence merely to pass a check. Share the relevant log/settings rather than transferring large matrices unnecessarily.
