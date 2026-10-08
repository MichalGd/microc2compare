# Resume acceptance after the compartment text-comparison failure

For a new installation from repository source, use [the GitHub installation guide](26-install-from-github.md). ZIP transfer sections here document earlier private deliveries; those artifacts are not included in Git.

Version **0.2.3** corrects the acceptance helper's byte-for-byte compartment comparison. It includes the existing TAD/modules, bootstrap correction and volcano plots. Compartment computation, statistical models, stored weights, environments and scientific completion-signature source files are unchanged from 0.2.2.

The user-reported 0.2.2 server outputs have identical coordinates, selected components and missing-value patterns for all six samples. Oriented-score differences are at most `1.55431e-15`. Several raw E1/E2/E3 vectors differ by a whole-vector sign; after sign alignment the largest difference is `2.88658e-15`. Thus the observed `cmp` failure reflects eigenvector sign freedom and floating-point rounding, not a demonstrated biological disagreement. Remaining acceptance checks still need to finish.

## What the corrected check requires

The dependency-free `scripts/compare_compartments.py` reads the two sample TSVs and their sibling `.eigenvalues.tsv`, `.phasing.json` and `.bedgraph` files. It checks:

- identical headers, row counts/order, genomic coordinates, selected components and other identity fields;
- identical missing-value patterns, with infinite/invalid values rejected;
- numerical agreement with fixed absolute tolerance `1e-10` and relative tolerance `1e-8`;
- at most one sign flip per raw component **per declared genomic view**, never independent per-bin flips or component swapping;
- oriented scores and browser values in the same direction, without sign alignment;
- unchanged view/eigenvalue identities, numerical eigenvalues, phasing status/component and correlation consistent with any raw sign flip.

The checker writes no scientific output and recalculates no quantities. Meaningful disagreements still exit nonzero. Its JSON stdout records tolerances, signs and maximum numerical errors. Insulation/loop acceptance comparisons retain their existing exact checks.

You can diagnose a single sample separately using the bootstrap Python:

```bash
"$m2c_bootstrap" "$m2c_src/scripts/compare_compartments.py" \
  --full "$m2c_accept/results/demo/synthetic_demo/run001/03_features/compartments/C1.tsv" \
  --selected "$m2c_accept/results/demo/synthetic_demo/compartments-only/03_features/compartments/C1.tsv"
```

Assign the variables in the installation steps below first. This command verifies the complete sample bundle, not just the displayed vector columns.

## 1. Transfer from Windows PowerShell

Run at the Windows `PS` prompt:

```powershell
$m2c_delivery = "D:\path\to\microc2compare\.local"
scp -P 22 "$m2c_delivery\microc2compare-0.2.3-agreement-8oct26.zip" "$m2c_delivery\microc2compare-0.2.3-agreement-8oct26.zip.sha256" "$m2c_delivery\install-microc2compare-0.2.3.py" "USERNAME@SERVER:/opt/conda_envs/microc2compare/incoming/"
```

## 2. Verify and install in server Bash

Retain the dedicated Bash/private cache setup from [candidate installation](20-candidate-upgrade.md#2-start-a-dedicated-linux-session). Restore these variables:

```bash
m2c_base="/opt/conda_envs/microc2compare"
m2c_work="$HOME/microc2compare-work"
m2c_bootstrap="/opt/miniconda/bin/python"
cd "$m2c_base/incoming"
sha256sum --check microc2compare-0.2.3-agreement-8oct26.zip.sha256
```

After checksum `OK`:

```bash
read -r m2c_sha m2c_archive_name < "$m2c_base/incoming/microc2compare-0.2.3-agreement-8oct26.zip.sha256"
"$m2c_bootstrap" "$m2c_base/incoming/install-microc2compare-0.2.3.py" \
  --archive "$m2c_base/incoming/microc2compare-0.2.3-agreement-8oct26.zip" \
  --sha256 "$m2c_sha" --base "$m2c_base" --release-id agreement-20261008
```

The installer verifies the manifest and refuses an existing release ID. After successful installation:

```bash
m2c_src="$m2c_base/releases/agreement-20261008/microc2compare-0.2.3"
bash "$m2c_src/scripts/microc2compare.sh" --help
"$m2c_bootstrap" "$m2c_src/scripts/compare_compartments.py" --help
```

This creates a separate source release, preserving earlier installations/archives and the shared command. No environment cloning, package installation or input generation is needed.

## 3. Resume the same acceptance fixture

```bash
m2c_site="$m2c_base/site/shared-design-json-fix-20261008.conf"
m2c_accept="$m2c_work/validation/modules-candidate-20261008"
set -o pipefail
bash "$m2c_src/scripts/validate_module_demo.sh" \
  --site "$m2c_site" --out "$m2c_accept" --resume \
  2>&1 | tee "$m2c_work/audit/modules-candidate-agreement-20261008.log"
m2c_accept_status=${PIPESTATUS[0]}
printf 'Candidate acceptance exit code: %s\n' "$m2c_accept_status"
```

The helper verifies original input hashes before execution. It refreshes preflight and reuses compatible completed science/report products. The release changes acceptance tooling only; matching scientific source digests are verified during local packaging. Ordinary signature checks still govern reuse if inputs/configuration/environments changed. Existing fixture/configurations are retained; missing selected configurations are created as the helper continues.

Success requires six compartment bundle checks, remaining selections/TAD variant, output validations and unchanged original matrices, ending with `Full/selected/variant acceptance passed` and exit code zero. Retain the log if any check fails. Do not edit completion records or bypass the numerical comparison. The demo remains descriptive; [volcano plots](22-loop-volcano-plots.md) require actual inferential model results.

## 4. Publish after acceptance and ordinary-account checks

Repeat acceptance from an ordinary researcher's own private workspace. After successful acceptance/report review, the maintainer creates the new versioned launcher:

```bash
"$m2c_bootstrap" "$m2c_src/scripts/publish_shared_command.py" --source "$m2c_src" --base "$m2c_base"
"$m2c_base/bin/microc2compare-0.2.3-agreement-20261008" --help
```

Then activate the accepted source explicitly:

```bash
"$m2c_bootstrap" "$m2c_src/scripts/publish_shared_command.py" --source "$m2c_src" --base "$m2c_base" --activate
hash -r
microc2compare --help
```

Keep previous versioned commands/releases. The existing `/usr/local/bin` symlink continues to work. Analysis/configuration/cache/scratch files stay in user-selected owned directories. Actual R/design calibration, biological pilots and GitHub publication remain separate pending work.
