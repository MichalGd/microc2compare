# Your first microc2compare analysis

microc2tracks converts paired sequencing reads into contact maps. microc2compare compares **existing maps from biological replicates** and writes compact tables, figures, reports and browser signals. You do not need to export a second balanced matrix or rerun alignment/balancing for an ordinary comparison.

## Understand your samples first

Suppose you have three control cultures and three treated cultures, all independently established. One control library was sequenced twice. That is six biological observations, not seven. The two sequencing sources belong to one library and require joint deduplication upstream. If the same culture instead yielded two independently prepared libraries, they remain one biological observation but can be deduplicated separately before consolidation.

For an unpaired experiment, give each actual culture a distinct `subject_id`, such as `culture_C1` and `culture_T1`. For matched samples from donor D1 before/after treatment, use the same actual donor ID in both conditions. Numbers like replicate `1` in two conditions do not establish pairing. [Inputs](02-inputs.md) explains the exact tables.

| Available biological observations | Eligible path under this workflow's policy |
|---|---|
| One per condition | Exploratory QC/features and descriptive effects; no biological differential P/FDR |
| Two independent units per condition | Descriptive, or explicitly preliminary if design and matching calibration pass |
| Two matched donor pairs | Descriptive; paired inference fails the required residual information floor |
| Three or more independent units per condition, or three matched pairs | Formal may be requested, subject to input/design/feature/calibration checks |

These are scoped eligibility policies, not guarantees of adequate power. Extra batch/covariate terms can make an otherwise eligible design unidentifiable or too small. Technical libraries, pools and downsampling repeats never increase biological sample size.

## Start with the synthetic demonstration

On an installed shared server, first follow [ordinary-user setup](16-shared-server-installation.md#3-start-as-an-ordinary-server-user). The maintainer follows [step-by-step server installation](15-server-installation-start.md) once to select compatible core and HiCRep prefixes. The [environment policy](01-installation.md) explains reuse and locks. Set the actual published release/profile paths below. Run from any working directory; the demo destination belongs to your account, outside the shared checkout:

```bash
m2c_src="/opt/conda_envs/microc2compare/releases/shared-design-json-fix-20261008/microc2compare-0.1.0"
m2c_site="/opt/conda_envs/microc2compare/site/shared-design-json-fix-20261008.conf"
source "$m2c_site"
export PYTHONDONTWRITEBYTECODE=1
umask 077
m2c_demo="$HOME/microc2compare-work/validation/first-demo-20261007"
"$CONDA_EXE" run --no-capture-output -p "$CORE_ENV_PREFIX" \
  python "$m2c_src/scripts/generate_demo.py" --out "$m2c_demo" \
  --core-prefix "$CORE_ENV_PREFIX" --raw-only
printf '\nsource "%s"\n' "$m2c_site" >> "$m2c_demo/project.conf"
```

The appended site profile selects the checked shared software without installing anything. The generated file sources the demo resource preset; reduce global budgets **and** enabled task reservations if the host allocation is smaller than four CPUs/16,000 decimal MB. `--raw-only` here disables cooltools features but retains weighted QC, HiCRep and reporting. It is a small first verification run, not a no-weight analysis and not formal inference. Choose a new empty generator destination if it already contains files.

```bash
bash "$m2c_src/scripts/bootstrap_environments.sh" -c "$m2c_demo/project.conf" --check-only
bash "$m2c_src/scripts/preflight_check.sh" -c "$m2c_demo/project.conf"
bash "$m2c_src/scripts/microc2compare.sh" -c "$m2c_demo/project.conf" --dry-run
bash "$m2c_src/scripts/microc2compare.sh" -c "$m2c_demo/project.conf"
```

The plan should show QC/HiCRep and reporting enabled, with pool/depth/features/differential modules disabled for this first demonstration. Expected outputs are under:

```text
$HOME/microc2compare-work/validation/first-demo-20261007/results/demo/synthetic_demo/run001/
  00_metadata/                  input/config/environment registrations
  01_qc/contacts/                six sample contact summaries
  01_qc/ps/                      raw and balanced contact-decay tables
  01_qc/hicrep/summary.tsv        fifteen biological-sample pairs
  08_reports/index.html          report with plots and downloads
  logs/                         native-tool/task diagnostics
  completion/                   validated completion records
```

Open `08_reports/index.html`. Check that all six samples appear, enabled modules have completed, and disabled modules are visibly skipped. Inspect a contact-decay plot and the concordance table. Undefined correlations remain missing values rather than zero. This checks file handling and software behavior on known synthetic data; it does not validate the biological method.

Check resume and output integrity:

```bash
bash "$m2c_src/scripts/microc2compare.sh" -c "$m2c_demo/project.conf" --resume
"$PYTHON_BOOTSTRAP" "$m2c_src/scripts/validate_outputs.py" \
  --config-json "$m2c_demo/results/demo/synthetic_demo/run001/00_metadata/config.json"
```

For the full descriptive feature demonstration, generate into a **different empty directory** without `--raw-only`, configure the checked prefixes, and repeat these steps. Full feature execution still requires the Linux/cooltools validation gates; the local verification record states what has actually run.

## Move to your upstream outputs

Use [the input helper and worked handoff](13-upstream-handoff.md) to inventory your real microc2tracks sample sheet/results and export draft manifests. Finish biological, library, deduplication, filtering and reference evidence. Start with explicit descriptive comparisons; passing a file scan does not authorize formal inference.

Your ordinary command remains:

```bash
bash "$m2c_src/scripts/microc2compare.sh" -c /your/project/config/project.conf --dry-run
bash "$m2c_src/scripts/microc2compare.sh" -c /your/project/config/project.conf
```

Use absolute server paths while learning. Quoting matters for paths with spaces. The `.conf` is trusted Bash; CSV/TSV rows are data and never shell commands. Results use `OUTDIR/PROJECT/ANALYSIS_SET_ID/RUN_ID`; change the run ID for a new historical run.

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
