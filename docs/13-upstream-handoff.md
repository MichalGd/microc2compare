# Preparing inputs from microc2tracks

For ordinary use, set `UPSTREAM_INPUTS` in your project configuration and run `microc2compare -c project.conf --prepare-inputs`; [the simple project guide](32-simple-project-commands.md) combines multiple upstreams without software paths or repeated project/reference/resolution arguments. The standalone `scripts/prepare_microc2compare_inputs.py` examples below are advanced alternatives for custom inventories before writing manifests by hand. It reads the upstream sample sheet, relates technical rows to biological units, suggests existing matrices, and lists associated outputs. It neither runs the upstream workflow nor alters its configuration, matrices, weights or contacts.

## Preview first

Run from the microc2compare checkout. Replace `/data/...` with your actual paths. A Python from the checked upstream/core environment normally already has `h5py` and NumPy:

```bash
python scripts/prepare_microc2compare_inputs.py \
  --samplesheet /data/upstream/config/samplesheet.csv \
  --upstream-results /data/upstream/results \
  --resolutions 10000,100000 --preview
```

The same command without `--preview` is also read-only. For a checked prefix without activation:

```bash
conda run --no-capture-output -p /actual/checked/core-prefix \
  python scripts/prepare_microc2compare_inputs.py \
  --samplesheet /data/upstream/config/samplesheet.csv \
  --upstream-results /data/upstream/results --preview --json
```

Human-readable output lists one suggested matrix per unit, constituent upstream rows, all resolution group names, requested-resolution integer count storage, weight validity/convergence, issues and related files. `--json` prints the full inventory. `needs_review` means files passed these limited inspection checks but biological/provenance review is still required. `blocked` means a missing, incompatible, ambiguous or uninspectable input was found. There is deliberately no automatic `ready_for_formal_analysis` status.

Exit code 0 means discovery completed, even when units are blocked; inspect each status. Invalid arguments/tables or an unsafe export return 2. A Python without `h5py`/NumPy can still discover file paths, but inspection is explicitly unavailable and units are blocked until inspected using a suitable Python.

## How matrices are selected

The authoritative biological key is reference, assay, condition and biological-replicate label. Reference aliases `mouse`/`GRCm39` map to `mm39`, and `human`/`GRCh38` map to `hg38`. This does not prove independence; subject/culture identity is reviewed separately.

When present, these existing files inform selection and membership checks:

```text
results/run_metadata/sample_manifest.tsv
results/run_metadata/technical_replicates.tsv
results/run_metadata/merge_manifest.tsv
results/final_report/microc2tracks_file_manifest.tsv
```

Without a usable manifest-listed matrix, the helper checks the established layout:

| Sources in the biological unit | Expected matrix |
|---|---|
| One upstream row | `results/{sample}/04_matrices/{sample}.norm.mcool` |
| Two or more upstream rows | `results/merged/{merge_group}/04_matrices/{merge_group}.norm.mcool` |

The default group name is `{condition}_{assay}_B{biological_replicate}_tech_merged`; explicit group names in the run manifest/sheet take precedence. Missing merged maps never fall back to individual technical matrices. A sheet containing only some members of a recorded merge is flagged. Conflicting group identities, duplicate matrix selections and mismatched run metadata are flagged as well. Manifest paths pointing to another run root are not selected; moved runs can use the expected layout.

Retained `.norm.mcool` files contain both raw pixels and stored balancing weights. Separate `.raw.mcool` files are not needed and may already have been removed upstream. The helper opens maps in read-only mode, reads Cooler headers, and scans requested bin-weight arrays in bounded chunks. It does not scan contact pixels, hash large matrices/FASTQs, compare all bin axes, establish count uniqueness, prove symmetry of square storage, or validate statistical designs. Normal preflight and upstream evidence still apply. See the [Cooler layout](https://cooler.readthedocs.io/en/latest/schema.html).

## Worked example: independent libraries and a single source

The [example sample sheet](../examples/upstream-handoff/samplesheet.csv) describes six independent cultures (three control, three treated). Culture C1 has two independently prepared libraries; the other cultures each have one upstream source. These example FASTQ paths are placeholders, not supplied data.

| Upstream rows | Actual biological unit | Suggested downstream sample/matrix |
|---|---|---|
| `C1_A`, `C1_B` | One culture, two independent libraries | `control_microc_B1_tech_merged` and its merged `.norm.mcool` |
| `C2` | Another independent control culture | `C2` and its individual `.norm.mcool` |
| `C3`, `T1`, `T2`, `T3` | Four additional independent cultures | One individual matrix per row |

Consequently seven FASTQ rows become six biological observations. Three is still the number of control biological replicates. Actual donors/cultures must be recorded as `subject_id`; the helper never derives those identities from condition/replicate labels.

For **lanes or resequencing of one library**, the membership should instead share the real `library_id`, use `source_kind=lane` or `resequencing`, and provide consistent joint-library deduplication evidence. Current upstream merging happens after per-row deduplication. If the lanes were independently deduplicated, reprocess them jointly upstream rather than filling in `joint_library` to bypass a check. Separate libraries can each be deduplicated independently and then consolidated within their biological unit. [Upstream merge implementation](https://github.com/MichalGd/microc2tracks/blob/main/scripts/merge_replicates.sh).

Helper-only optional sheet fields `subject_id`, `batch`, `library_id`, `source_kind`, `dedup_group_id`, `dedup_scope` and `source_checksum` can carry explicitly reviewed information forward. Their absence stays visible. The helper reads the original sheet without rewriting it; you can instead complete the exported drafts.

## Export editable drafts

### Older manifests and deliberately selecting individual maps

Version 0.2.7 also reads older sample/merge manifests without `reference_id` and the legacy headerless `technical_replicates.tsv` layout: `merge_group`, `sample`, valid-pairs path. It retains missing reference evidence; the confirmed reference must come from a sheet/run record or explicit `--legacy-reference-id`. The option never overrides a known conflicting reference. Malformed records still fail; upstream files are never rewritten.

An older run may have pooled all samples even though your reviewed sheet now describes separate biological units. To deliberately select **one existing unmerged sample map per declared biological unit**, add:

```bash
python scripts/prepare_microc2compare_inputs.py \
  --samplesheet /data/reviewed/samplesheet.csv \
  --upstream-results /data/upstream/results \
  --legacy-reference-id mm39 --individual-matrices --preview
```

Each sheet biological unit must have exactly one source. Multiple technical inputs still require upstream consolidation; this flag does not manufacture biological replication or split a pooled matrix. The same flag applies to `--write-drafts`. Historical biological/technical labels and pooled-group differences stay visible for review, and the old sample/technical records are copied to the inventory and provenance sidecars. Missing sample-manifest records are explicitly reported. Actual assay/condition/reference conflicts, missing/invalid maps, duplicate selections and shared subject/library errors still block.

The default behavior without this flag retains strict historical membership checks, including partial-merge rejection. Use it only after checking that the individual maps represent the units in your reviewed sheet. Check count uniqueness, experimental identities and processing evidence before clearing draft fields.

Choose a **new directory outside upstream results**:

```bash
python scripts/prepare_microc2compare_inputs.py \
  --samplesheet /data/upstream/config/samplesheet.csv \
  --upstream-results /data/upstream/results \
  --project experiment --reference-id mm39 \
  --write-drafts /data/experiment/config
```

Do this before creating `/data/experiment/config`; existing destinations are refused, even if empty. Exporting remains a preparation step; blocked units can be recorded for diagnosis. Files are:

```text
config/
  READ_BEFORE_RUNNING.md   issues and mandatory review
  input_inventory.json    complete machine-readable discovery
  input_inventory.tsv     sortable unit/status summary
  related_outputs.tsv     pairs, QC, .hic and other discovered files
  samples.csv             one suggested biological-unit row, comma-separated
  technical_membership.tsv
  comparisons.tsv         header only; choose contrasts explicitly
  references.json         reference-registration skeleton
  provenance/*.json       actual stored weight evidence plus unresolved fields
```

Paths to imported matrices are absolute. Provenance paths are relative to `samples.csv`. FASTQ strings in the membership table are copied verbatim as upstream evidence; resolve relative FASTQ paths against the original upstream setup before using them to locate or checksum sources.

Complete the following before running:

1. Confirm that grouping reflects actual biological units; fill `subject_id` and needed batch/covariates. Matched donors reuse a subject ID across conditions and require an explicit paired comparison.
2. Complete real library/source identities, source kinds, deduplication scopes/groups and original-source checksums. Never substitute the matrix checksum for a production source checksum.
3. Review each provenance JSON: declare actual count units, filtering, processing version and deduplication. Stored convergence is copied as true/false/null; missing evidence is never changed to true. Recipe parameters are copied only when present in weight attributes. Verify the complete compatible upstream recipes; an empty or partial recipe is not acceptance evidence. Resolve null filtering fields and the `draft` marker after review.
4. Register assembly-matched chromosome sizes, analysis views, their SHA256 values, and the phasing track if compartments are enabled. Upstream reference TSVs are listed where found, but they are not automatically turned into an accepted downstream registry. See [input contracts](02-inputs.md).
5. Add explicit rows to `comparisons.tsv`; start with `inference=descriptive`. The helper deliberately does not choose a control, subject pairing, design or formal status.

For a mixed-reference upstream run, preview everything, then export one reference at a time using `--reference-id`. Export refuses mixed reference identities. A legacy sheet with missing reference fields needs actual run metadata or your explicit `--legacy-reference-id mm39`; there is no silent mouse default.

After review, copy `config/config_template.conf` to the new configuration directory as `project.conf`. Set `PROJECT`, `REFERENCE_ID`, `SAMPLESHEET="samples.csv"`, `TECHNICAL_MEMBERSHIP="technical_membership.tsv"`, `COMPARISONS="comparisons.tsv"`, `REFERENCE_REGISTRY="references.json"`, environment prefixes and output/scratch paths. Then:

On a shared server, this configuration directory belongs in your user/project workspace, outside the shared source. The profile-bound 0.2.5 command selects software internally; no manual software source line is needed. Do not create environments or edit the shared template. See [ordinary-user setup](16-shared-server-installation.md#3-start-as-an-ordinary-server-user).

```bash
microc2compare -c /data/experiment/config/project.conf --dry-run
# Review the plan and validation status before starting computation.
microc2compare -c /data/experiment/config/project.conf --until-stage qc
```

## Compatibility choices

Keep the upstream-style Bash runner, `.conf` configuration, `OUTDIR`, module/thread/concurrency switches, CSV/TSV tables and numbered guides. Retain the existing stage folders and testable Python package. The [ATAC workflow](https://github.com/MichalGd/ATACseq2tracks) provides useful planning, recovery and replicate-audit patterns; the [RNA workflow](https://github.com/MichalGd/rna_ends2tracks) provides useful metadata, report and documentation-index patterns.

Unlike their FASTQ sheets, the downstream sheet contains consolidated biological observations and needs separate technical membership/reference/contrast evidence. Comparisons require common references/grids and explicit discovery/testing separation. Pools remain descriptive discovery inputs, not model replicates. Unlike newer literal-only ATAC/RNA configuration parsers, microc2compare currently sources a trusted Bash `.conf`, including explicit resource presets. Do not source the upstream configuration to discover outputs; the helper never executes it.

## Common problems

| Observation | Action |
|---|---|
| Merged map missing but individual maps exist | Complete/review upstream consolidation; never count the technical rows as independent replicates |
| Unknown reference | Supply actual run metadata or explicitly declare the confirmed legacy reference |
| Missing resolution or failed balancing | Correct upstream output; do not rebuild or rebalance imported samples downstream |
| Convergence unknown | Retrieve actual upstream evidence or review the explicitly documented legacy exception; do not invent convergence |
| Partial/mismatched merge membership | Use the corresponding full sheet/run or reconstruct the correct biological-unit map upstream |
| Multiple biological identities share a group/library | Correct actual identity/grouping; do not merely rename rows to manufacture replication |
| Export directory exists | Choose a new destination and compare drafts with the existing project manually |
| Matrix inventory looks fine but preflight fails | Inspect its specific evidence/design/reference/environment failure; discovery is deliberately less comprehensive |
