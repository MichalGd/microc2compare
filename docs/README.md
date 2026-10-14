# Documentation index

For the current 0.4.0 source candidate, start with the [general-purpose upgrade guide](35-general-purpose-upgrade.md). It covers new projects, CPU/RAM admission, fine grids, stripes, browser exports and migration from an accepted release. Project inputs and references come from the current manifests; no earlier pilot cohort is required. Numbered filenames remain stable, and earlier guides retain useful setup and release history. Original material in `assumptions/` and `prompts/` records decisions rather than current run instructions.

## First use

| Question | Guide |
|---|---|
| How do I configure a new project and install the 0.4.0 candidate alongside an accepted release? | [General-purpose upgrade](35-general-purpose-upgrade.md) |
| How do I prepare/run/check a project through the shared command? | [Simple project commands](32-simple-project-commands.md) |
| How do I install the source from GitHub on a shared server? | [Repository installation](26-install-from-github.md) |
| How do I run the shared command without a manual software-profile line? | [Automatic software settings](25-automatic-software-profile.md) |
| What is a matrix, loop or biological replicate? How do I run the demonstration? | [First analysis and glossary](14-first-analysis.md) |
| How do I select the right microc2tracks outputs? | [Upstream handoff](13-upstream-handoff.md) |
| How do I install one shared clone without changing home/upstream software? | [Step-by-step server installation](15-server-installation-start.md) |
| How do I use the shared installation from my own account? | [Shared installation: user setup](16-shared-server-installation.md#3-start-as-an-ordinary-server-user) |
| What are the environment reuse/lock policies? | [Linux installation](01-installation.md) |
| What belongs in sample, membership and reference files? | [Inputs](02-inputs.md) and [current observation identity](35-general-purpose-upgrade.md#observation-identity-and-migration) |
| Which settings and analyses should I choose? | [Configuration](03-configuration.md), [module catalogue](17-analysis-modules.md) and [current upgrade](35-general-purpose-upgrade.md) |

## Current analyses and exports

| Question | Guide |
|---|---|
| How do the generic eight-grid and short-range presets work, and what can target 4 kb genomic separation? | [CPU loop profiles, caller geometry and input capability](cpu_loop_profiles.md) |
| How do I prepare authentic missing grids from retained exact-position pairs? | [Explicit pairs preparation](cpu_loop_profiles.md#explicit-exact-position-pairs-preparation) |
| How does independent Chromosight discovery differ from fixed-candidate scoring? | [CPU loop methods](cpu_loop_profiles.md) |
| How do I detect stripes and interpret native StripePy descriptors? | [Stripe detection](34-stripe-detection.md) |
| How do I export saved individual/condition tracks to UCSC or optional HiGlass formats? | [Unified browser export](browser_export.md) |
| How do I measure TAD growth, merging or disappearance? | [TAD analysis and event outputs](18-tad-analysis.md) |
| How do I annotate both loop anchors? | [Gene/cCRE annotation](29-loop-anchor-annotation.md) |
| How do I select autosomes and optional X/Y? | [Chromosome scope](30-chromosome-scope.md) and [common chromosome alignment](27-common-chromosomes.md) |
| How do I interpret differential loop volcano plots? | [Volcano plots](22-loop-volcano-plots.md) |

The new loop presets adapt published methods and keep scientific admission project-specific. They do not guarantee recovery at 4 kb or reproduce the historical Hsieh study. Multiple callers from one observation do not add biological replication. Stripe support and condition browser summaries remain descriptive.

## Routine operation

| Question | Guide |
|---|---|
| How do I preview, run, monitor or resume? | [Operation](06-operation.md) and [current checkpoint/migration behavior](35-general-purpose-upgrade.md) |
| Where are results, and what do they mean? | [Outputs and worked interpretation](09-outputs.md), [browser export](browser_export.md) and [stripe outputs](34-stripe-detection.md) |
| How do I diagnose a failed check? | [Troubleshooting](10-troubleshooting.md) |
| How are CPU/RAM and matrix-reader limits handled? | [Current resource admission](35-general-purpose-upgrade.md) and [server resources](07-server-resources.md) |

## Methods, statistics and validation

| Question | Guide |
|---|---|
| How are counts, weights, loops, insulation and aggregates calculated? | [Methods](04-methods.md) and [CPU caller/pairs details](cpu_loop_profiles.md) |
| When is differential inference allowed? How are pairing and batch handled? | [Statistics](05-statistics.md) and [donor/culture limitations](35-general-purpose-upgrade.md#observation-identity-and-migration) |
| Which 0.4.0 checks are local, and which Linux/biological checks remain pending? | [Current validation boundaries](35-general-purpose-upgrade.md#validation-boundaries) |
| Where is earlier validation evidence retained? | [Validation guide](08-validation.md) and [verification record](../validation/local-verification.md) |

Completed local checks cover synthetic new projects, numerical agreement, input preservation, recovery and actual caller APIs, including isolated StripePy execution. The final 0.4.0 full-suite tally is pending. Linux CI/converter execution, shared-server resource acceptance, real-data sensitivity, installed-browser rendering and formal biological calibration remain separate checks. Earlier accepted releases do not establish acceptance of this candidate.

## Administration and development

| Question | Guide |
|---|---|
| How do I deploy and accept 0.4.0 without replacing running code? | [General-purpose upgrade and migration](35-general-purpose-upgrade.md) |
| How do I make the workflow available to every server account? | [Shared server installation and access tests](16-shared-server-installation.md) |
| Can I reuse an R environment, and what is required for formal inference? | [R environment and statistical acceptance](19-r-environment-and-inference.md) |
| How do I prepare a reviewed source archive and GitHub release? | [Publication](11-publication.md) |
| What works on Windows, and which checks require Linux? | [Windows development](12-windows-development.md) and [current validation boundaries](35-general-purpose-upgrade.md#validation-boundaries) |
| How do I contribute or run checks? | [Contributing](../CONTRIBUTING.md) |

Historical release-specific procedures remain available: [0.3.2 convenience upgrade](33-simple-project-upgrade.md), [fine-loop upgrade](31-fine-loop-upgrade.md), [2/4-kb loops and required saddles](28-fine-resolution-loops.md), [module/TAD candidate upgrade](20-candidate-upgrade.md), [bootstrap correction](21-bootstrap-state-correction.md), [compartment agreement correction](23-compartment-agreement.md) and [R loop correction](24-r-loop-compatibility.md). Use the current guide for 0.4.0 installation and validation status; legacy full practical/Hsieh-inspired profiles retain their required compartment/saddle/report behavior.

The [README](../README.md) contains the workflow diagram and shortest command sequence. Formal inference remains subject to the documented design and calibration gates.
