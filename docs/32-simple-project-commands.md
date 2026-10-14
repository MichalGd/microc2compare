# Simple project commands (0.3.2)

Use your normal Linux Bash login and the installed `microc2compare` command. No Conda activation, software-profile `source`, versioned source path, Python path, cache exports or separate clean Bash session is needed for these commands. The shared launcher chooses software; the runner sets private permissions, removes inherited Python/R library overrides and creates runtime caches under `SCRATCH_ROOT/.runtime`.

## 1. Create a small configuration

```bash
microc2compare --init "$HOME/my_microc_project/config/project.conf"
```

This creates a compact `project.conf`, an empty `upstream_inputs.tsv`, an empty `comparisons.tsv` and a reference descriptor. It refuses existing files. Edit the experiment-specific values: project/run identifiers, genome, input locations, desired outputs/scratch and your CPU/RAM budget. `OUTDIR="../results"` and `SCRATCH_ROOT="../scratch"` are relative to the configuration directory; both remain inside your project in this example. They can instead name any appropriate directory you own. No results go into the shared software installation.

The generated `use_analysis_profile practical_2_4kb` line loads versioned defaults for 2/4-kb Mustache loops, 10-kb insulation/TAD/ADA, 100-kb compartments, APA and obligatory saddle plots in the full report. It is supplied automatically, with no path to fill in. Put optional parameter overrides **below** it. Software requirements and scientific eligibility gates still apply. Use `--init ... --preset existing` for the older 10-kb profile. Existing full configurations remain supported; initialization is unnecessary for an already prepared project.

To try the software before adding real inputs, `microc2compare -c project.conf --demo` creates a separate six-sample synthetic QC fixture and prints its configuration path. Run the printed command and then `--validate` with that fixture's configuration. This uses the shared core/HiCRep software automatically and never replaces your real configuration or data. It is a software demonstration, not biological validation.

## 2. Prepare inputs only for a new upstream handoff

If you already have completed `samples.csv`, technical membership and provenance, point the configuration to them and skip this step. Otherwise edit `upstream_inputs.tsv` with one row per upstream dataset. Its required columns are `samplesheet` and `upstream_results`; `individual_matrices` is an optional explicit `true`/`false` choice. Paths are relative to this table unless absolute.

| samplesheet | upstream_results | individual_matrices |
|---|---|---|
| /your/project/mESCs.samples.csv | /your/project/mESCs-regrouped | false |
| /your/project/basal.samples.csv | /your/upstream/basal | true |

`false` selects matrices for the biological units described in the upstream sheet, including already merged technical sources. `true` explicitly selects individual matrices for the sheet's biological units, retaining conflicting historical metadata for review. It never changes technical rows into biological replicates or silently merges maps. For this project's mESCs, rep1/rep2 are technical members of each declared cell-line unit; basal MC1–MC4 remain four biological observations.

```bash
microc2compare -c "$HOME/my_microc_project/config/project.conf" --prepare-inputs
```

Project, reference identity and requested resolutions come from the configuration. The command inventories all listed upstreams, combines their sample and technical membership tables and preserves separate evidence inventories. Draft provenance remains explicitly incomplete where experimental facts are absent. Read the printed `READ_BEFORE_RUNNING.md` and complete genuine library/deduplication/filtering evidence; passing discovery is not preflight acceptance. Missing historical reference labels use your explicitly declared `REFERENCE_ID`, while conflicts and evidence gaps remain visible.

Preparation prints a **new sibling configuration** pointing to these drafts. Use that file for the next steps. The original configuration and upstream outputs are preserved. No source checksum scanning, balancing or large matrix copying is triggered by input discovery.

## 3. Register reference and optional annotation assets once

Reuse a checked reference registry announced by the maintainer when appropriate. Users sharing the same immutable assets can share its read-only registry. They do not need to register or copy it again per analysis. Otherwise set `REFERENCE_REGISTRY` to a JSON/YAML descriptor with the genuine `reference_id` and `chromsizes` path. Include existing `views`, `blacklist` and `phasing_track` paths when applicable. A common phasing track is a signal used to orient compartment vectors consistently; its presence and usability are required for compartments and saddles.

```json
{
  "reference_id": "mm39",
  "chromsizes": "/actual/reference/mm39.chrom.sizes",
  "phasing_track": "/actual/reference/mm39.100kb.phasing.bedgraph"
}
```

Replace these two asset paths with files you actually have. Do not enter SHA256 strings by hand. If annotation is wanted, set `ANNOTATION_REGISTRY` to a descriptor containing the actual assembly, GTF path/release and optional cCRE path/source/conversion history. Release and conversion information are scientific facts, not invented defaults. Then run:

```bash
microc2compare -c /your/project/config/project.conf --register-resources
```

The command calculates missing SHA256 hashes, resolves paths, validates reference coordinates and writes new registries plus a sibling configuration. If views are absent, it generates whole-chromosome views from the size file; the analysis scope still defaults to autosomes. Existing valid hashes are checked; mismatched hashes fail rather than silently accepting changed assets. Original assets are never modified or copied. No phasing or annotation signal is invented, and this operation does not certify annotation semantics or statistical calibration. Full preflight remains mandatory. Set `RUN_ANNOTATION=true` only when the chosen annotation inputs are ready.

## 4. Audit or explicitly prepare requested resolutions

```bash
microc2compare -c /your/project/config/project.conf --audit-resolutions
```

Resolutions, samples, registry overrides and chromosome policy come from the configuration. `availability.tsv` and `audit.json` report missing grids, weight convergence, positive usable weights and metadata identities. The audit does not scan all contacts or approve coverage/reproducibility. To inspect an additional grid, add its actual scientific request in the configuration, or use the advanced standalone helper described in [fine-resolution analysis](28-fine-resolution-loops.md).

If a requested grid is missing and an exact finer raw grid exists:

```bash
microc2compare -c /your/project/config/project.conf --prepare-resolutions
```

This **explicit** command reuses existing grids (including registered derived grids), creates only missing integer coarsenings in a new project directory and calculates weights only for those new matrices. It verifies raw count conservation and unchanged original file hashes. It cannot reconstruct 2-kb data from a 10-kb map: that requires upstream contact-pair processing. It honors the configured CPU reservation, observed allocation and balancing memory admission budget; one matrix is prepared at a time. The reservation is not an operating-system RSS limit. A failure leaves incomplete products separate from accepted inputs and prints no ready configuration. This operation may take hours and is never automatic during ordinary analysis.

Successful preparation prints a new configuration containing the generated sample and matrix registry pointers. Users do not have to assemble those pointers themselves. Prepared-resolution balancing uses the actual stored upstream recipe and requires convergence; downstream analysis still never rebalances imported matrices.

## 5. Run and verify

```bash
microc2compare -c /your/project/config/project.conf --dry-run
microc2compare -c /your/project/config/project.conf --analysis qc,hicrep
microc2compare -c /your/project/config/project.conf --validate
```

`--validate` locates the run from `OUTDIR/PROJECT/ANALYSIS_SET_ID/RUN_ID`, checks the modules recorded for that run and compares current configuration/source/input identities and required output hashes. Missing or stale tasks give a nonzero exit code. Changing relevant parameters after a run can therefore invalidate it. Validate old runs with their original versioned launcher and configuration. Validation does not mean the biological method was calibrated.

For the full analysis use a new `RUN_ID`, check your references/fine-grid pilot and run without `--analysis`. For a selected module, choose another run ID or use reviewed `--resume` behavior. See [modules](17-analysis-modules.md), [outputs](09-outputs.md) and [statistics](05-statistics.md). Saddle plots remain obligatory in full fine-loop profiles and appear in the final HTML report; standalone QC/loops/other selections remain possible.

To change parameters in a separate history without overwriting the original configuration:

```bash
microc2compare -c /your/project/config/project.conf --variant tad_threshold_05 \
  --set INSULATION_BOUNDARY_THRESHOLD=0.5
```

The command prints a new sibling configuration with the new run ID and literal overrides. Dry-run that file, then run `--analysis tads` or your other chosen modules. Relative input/output paths retain their meaning. Repeated `--set` adds parameter choices; no separate Python/source path is needed. A variant does not run the analysis or automatically establish a statistically valid setting.

By default preparation/audit products go under `OUTDIR/PROJECT/preparation/ACTION-TIMESTAMP-ID`. Each operation uses a new directory. Optional `--out /your/chosen/new-directory` replaces that destination; it never changes the analysis `OUTDIR`. These actions are separate from execution options (`--resume`, `--dry-run`, `--analysis`, stage selection). Keep the printed configurations and inventories for reproducibility; delete disposable preparation outputs according to your own project's retention policy.

For jobs that must survive SSH disconnection, `nohup` remains an explicit execution choice:

```bash
nohup microc2compare -c /your/project/config/project.conf \
  > /your/project/analysis.log 2>&1 < /dev/null &
```

All paths above are experiment-specific. Fixed software/environment paths and cache setup are internal to the installed command. Installation, environment creation and calibration acceptance remain maintainer/method-validation procedures, not repeated ordinary-user setup.
