# Documentation index

Start with the first-analysis guide. The numbered filenames remain stable; their order below follows the questions a new researcher is likely to ask. Original planning material in `assumptions/` and `prompts/` records decisions and is not a replacement for current run instructions.

## First use

| Question | Guide |
|---|---|
| How do I install the source from GitHub on a shared server? | [Repository installation](26-install-from-github.md) |
| How do I run the shared command without a manual software-profile line? | [Automatic software settings](25-automatic-software-profile.md) |
| What is a matrix, a loop or a biological replicate? How do I run the demonstration? | [First analysis and glossary](14-first-analysis.md) |
| How do I select the right microc2tracks outputs? | [Upstream handoff](13-upstream-handoff.md) |
| How do I install one shared clone without changing home/upstream software? | [Step-by-step server installation](15-server-installation-start.md) |
| How do I use the shared installation from my own account? | [Shared installation: user setup](16-shared-server-installation.md#3-start-as-an-ordinary-server-user) |
| What are the environment reuse/lock policies? | [Linux installation](01-installation.md) |
| What do I put in the sample, membership and reference files? | [Inputs](02-inputs.md) |
| Which settings should I change? | [Configuration](03-configuration.md) |
| Which analyses exist, what inputs do they need, and how do I run one? | [Module catalogue and parameter variants](17-analysis-modules.md) |
| How do I measure TAD growth, merging or disappearance? | [TAD analysis and event outputs](18-tad-analysis.md) |
| How do I interpret/export differential loop volcano plots and install the amended candidate? | [Volcano plots and 0.2.2 installation](22-loop-volcano-plots.md) |

## Routine operation

| Question | Guide |
|---|---|
| How do I preview, run, monitor or resume? | [Operation](06-operation.md) |
| Where are my results, and what do they mean? | [Outputs and worked interpretation](09-outputs.md) |
| How do I diagnose a failed check? | [Troubleshooting](10-troubleshooting.md) |

## Methods and statistics

| Question | Guide |
|---|---|
| How are counts, weights, loops, insulation and aggregates calculated? | [Methods](04-methods.md) |
| When is differential inference allowed? How are pairing and batch handled? | [Statistics](05-statistics.md) |
| What has actually been validated? | [Validation guide](08-validation.md) and [verification record](../validation/local-verification.md) |

## Administration and development

| Question | Guide |
|---|---|
| How do I inventory and tune the 500-GB/70-core server? | [Server resources](07-server-resources.md) |
| How do I make the workflow available to every server account? | [Shared server installation and access tests](16-shared-server-installation.md) |
| How do I install/accept the new module/TAD release beside the working installation? | [Candidate upgrade](20-candidate-upgrade.md) |
| How do I resume after the 0.2.0 bootstrap missing-Cooler failure? | [Bootstrap state correction](21-bootstrap-state-correction.md) |
| Why do full/selected compartment TSVs differ, and how do I resume? | [Numerical compartment agreement correction](23-compartment-agreement.md) |
| Can I reuse an R environment, and what is required for formal inference? | [R environment and statistical acceptance](19-r-environment-and-inference.md) |
| Why does R loop validation fail with `vapply ... result is length 0`? | [R loop correction and 0.2.4 acceptance](24-r-loop-compatibility.md) |
| How do I prepare a reviewed source archive and later GitHub release? | [Publication](11-publication.md) |
| What works on Windows, and which checks still require Linux? | [Windows development](12-windows-development.md) |
| How do I contribute or run checks? | [Contributing](../CONTRIBUTING.md) |

The [README](../README.md) contains the workflow diagram and shortest command sequence. Formal inference remains subject to the documented validation gates.
