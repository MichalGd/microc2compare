# Publication and release preparation

The existing repository is `MichalGd/microc2compare`. Prepare upgrades in a separate checkout of its current history. A local source commit and an authenticated remote push are separate steps; verify both rather than treating preparation as publication. Publish only reviewed source, tests, examples, documentation, schemas, environment definitions and nonprivate validation records.

## Public archive

```bash
python3 scripts/package_release.py --out /owned/releases/microc2compare-0.4.0.zip
```

The archive is an explicit allowlist of source, tests, public examples, documentation, schemas, environment definitions and nonprivate validation records. It excludes original assumptions/prompts, environments, private configurations/manifests, contact matrices, generated results/scratch, caches, PDF files and secrets. Original workspace documents are preserved locally. Review the archive listing and every included file before publication; an allowlist cannot prove an arbitrary future file contains no sensitive content.

## Repository publication

1. Review source/scientific evidence and the pending release gates; avoid a production-ready claim.
2. Confirm visibility and licensing. A private repository can be prepared while licensing is undecided; public distribution requires the maintainer's explicit choice. [LICENSE.md](../LICENSE.md) records the current decision. No upstream license/owner is inferred.
3. Review private data, absolute paths, hashes, credentials and redistribution rights; publish only the explicit archive contents.
4. Initialize an owned Git checkout, inspect `git status`/diff, stage the public allowlist, and review the commit. Do not use `git add .` on a private analysis workspace without inspection.
5. Push the reviewed commit on an upgrade branch in the existing repository using authenticated maintainer access. Open a pull request so Linux validation results are visible. Refuse to overwrite another repository or force-push. Preserve the existing repository visibility and license.
6. Verify the URL, visibility, `main` commit and file set. Inspect actual Linux/R CI results; source publication does not certify hosted CI or statistical calibration. Tag only reviewed commits with accurate candidate status, and complete applicable calibration/real-data gates before production endorsement.

From an authenticated GitHub CLI account with repository write permissions:

```bash
git push -u origin upgrade/general-purpose-0.4.0
gh pr create --repo MichalGd/microc2compare --base main \
  --head upgrade/general-purpose-0.4.0 --title "General-purpose CPU and browser upgrade" \
  --body-file /owned/reviewed-upgrade-description.md
```

The [official CLI guide](https://cli.github.com/manual/gh_pr_create) documents pull-request creation. Authentication belongs in normal GitHub tooling, never source files or project configurations. The source release remains a validation candidate until the required Linux/target-server gates pass.

New installations follow [the repository installation guide](26-install-from-github.md). Historical ZIP transfer guides describe private deliveries that are not GitHub assets. Preserve their immutable archives/checksums; repository documentation can be updated without replacing those deliveries.

Dependencies retain their own licenses. No microc2tracks source was copied. Candidate public-data registrations require exact accessions and redistribution checks before including any data. Cite the scientific tools independently of this project.

The local release package is a reviewable artifact, not a legal license grant, deployment or public release. `CITATION.cff` gives the software title/version and a collective contributor entry; replace it with approved authors/identifiers when ownership and authorship are decided.
