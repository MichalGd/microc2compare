# Publication and release preparation

The maintainer has requested creation of `MichalGd/microc2compare`. A local source commit and actual remote creation/push are separate steps; verify both rather than treating preparation as publication. Publish only reviewed source, tests, examples, documentation, schemas, environment definitions and nonprivate validation records.

## Public archive

```bash
python3 scripts/package_release.py --out /owned/releases/microc2compare-0.2.5.zip
```

The archive is an explicit allowlist of source, tests, public examples, documentation, schemas, environment definitions and nonprivate validation records. It excludes original assumptions/prompts, environments, private configurations/manifests, contact matrices, generated results/scratch, caches, PDF files and secrets. Original workspace documents are preserved locally. Review the archive listing and every included file before publication; an allowlist cannot prove an arbitrary future file contains no sensitive content.

## Repository publication

1. Review source/scientific evidence and the pending release gates; avoid a production-ready claim.
2. Confirm visibility and licensing. A private repository can be prepared while licensing is undecided; public distribution requires the maintainer's explicit choice. [LICENSE.md](../LICENSE.md) records the current decision. No upstream license/owner is inferred.
3. Review private data, absolute paths, hashes, credentials and redistribution rights; publish only the explicit archive contents.
4. Initialize an owned Git checkout, inspect `git status`/diff, stage the public allowlist, and review the commit. Do not use `git add .` on a private analysis workspace without inspection.
5. Create the authorized GitHub repository using authenticated maintainer access; add the chosen remote and push the reviewed commit. Refuse to overwrite another repository or force-push. Do not initialize a competing GitHub README/license commit when pushing the prepared history.
6. Verify the URL, visibility, `main` commit and file set. Inspect actual Linux/R CI results; source publication does not certify hosted CI or statistical calibration. Tag only reviewed commits with accurate candidate status, and complete applicable calibration/real-data gates before production endorsement.

The private initial publication command, from an authenticated GitHub CLI account with ownership permissions, is:

```bash
gh repo create MichalGd/microc2compare --private \
  --description "Bash-based replicate-aware comparison of existing Micro-C/Hi-C matrices" \
  --source /path/to/reviewed/git-checkout --remote origin --push
```

Use `--public` only for an approved public/licensed version. The [official CLI guide](https://cli.github.com/manual/gh_repo_create) documents the options. Authentication belongs in normal GitHub tooling, never source files or project configurations.

New installations follow [the repository installation guide](26-install-from-github.md). Historical ZIP transfer guides describe private deliveries that are not GitHub assets. Preserve their immutable archives/checksums; repository documentation can be updated without replacing those deliveries.

Dependencies retain their own licenses. No microc2tracks source was copied. Candidate public-data registrations require exact accessions and redistribution checks before including any data. Cite the scientific tools independently of this project.

The local release package is a reviewable artifact, not a legal license grant, deployment or public release. `CITATION.cff` gives the software title/version and a collective contributor entry; replace it with approved authors/identifiers when ownership and authorship are decided.
