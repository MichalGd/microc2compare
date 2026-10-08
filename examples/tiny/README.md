# Tiny original demonstration

`scripts/generate_demo.py` creates six redistributable synthetic biological samples in an empty destination: three control and three treated cultures, a 6-Mb chromosome, 10/100-kb groups, decay/compartment patterns, a focal contact change, immutable reference assets and complete descriptive manifests. Known unit test weights are stored upstream in each fixture map; the workflow does not calculate them.

Generate with a checked core prefix on Linux and set the checked HiCRep prefix in the generated `.conf`. Use the demo resource preset or lower both budgets and task reservations on a smaller allocation. See [README](../../README.md) for the complete command and [Windows development](../../docs/12-windows-development.md) for the actual local raw-only commands.

Default pooling, thinning and new-map balancing remain off. Full mode exercises actual cooltools descriptive features; `--raw-only` disables those and differential/aggregate modules but retains weighted QC, HiCRep and reporting. No formal calibration claim follows from this synthetic example. The generator refuses to overwrite a nonempty destination; generated maps/results are ignored and excluded from the public archive.
