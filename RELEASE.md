# Publishing ctgimme 0.0.12

Version 0.0.12 is the CRAN resubmission release. It renames both the package
and public functions to `ctgimme`, addresses CRAN's example and console-output
comments, and retains the multisubject subgroup likelihood introduced in the
earlier 0.0.11 package version.

The existing local Git tag `v0.0.11` names the prior `ctsgimme` package
version. Do not move or overwrite it. Publish it together with the historical
`pre-package-0.0.6` tag if those tags are still absent from the remote.

## Accepted archive

The final checked source archive is `dist/ctgimme_0.0.12.tar.gz` (51,021
bytes), SHA-256
`485E698808387C8FB15E30EC9B464878BCF88B82BEF710ED2D756CE955AF1ABA`.
On Windows 11 with R 4.6.1, `R CMD check --as-cran --run-donttest
--timings` completed with 0 errors, 0 warnings, and two local notes: the
expected new-submission note and skipped HTML validation because the local
machine lacked an HTML Tidy executable. The package example elapsed time was
0.68 seconds, and all 190 test expectations passed.
The retained log, timing file, and installed snapshot for this exact archive
are under `dist/cran-final-0.0.12-h/ctgimme.Rcheck/`; older candidate check
directories are not release evidence.

## Local release checks

Run these commands from the repository root with current R release or
R-patched, a working LaTeX installation, and Pandoc on `PATH`:

```powershell
Rscript --vanilla -e "roxygen2::roxygenise('.')"
Rscript --vanilla -e "devtools::test(stop_on_failure = TRUE)"
Rscript --vanilla -e "spelling::spell_check_package('.')"
Rscript --vanilla -e "urlchecker::url_check('.')"

New-Item -ItemType Directory -Force dist | Out-Null
Push-Location dist
R CMD build ..
R CMD check --as-cran --run-donttest --timings ctgimme_0.0.12.tar.gz
Pop-Location
```

Inspect `dist/ctgimme.Rcheck/00check.log` and
`dist/ctgimme.Rcheck/ctgimme-Ex.timings`. The main `ctgimme()` example must
remain below five seconds; it is intentionally unwrapped and writes only to a
temporary directory.

Install and smoke-test the exact archive in a new library:

```powershell
$releaseLibrary = Join-Path $env:TEMP ("ctgimme-0.0.12-" + [guid]::NewGuid())
New-Item -ItemType Directory -Force $releaseLibrary | Out-Null
R CMD INSTALL "--library=$releaseLibrary" dist/ctgimme_0.0.12.tar.gz
if ($LASTEXITCODE -ne 0) { throw "R CMD INSTALL failed." }
$effectiveLibraries = Rscript --vanilla -e 'cat(paste(.libPaths(), collapse = .Platform$path.sep))'
if ($LASTEXITCODE -ne 0) { throw "Reading the R library path failed." }
$priorLibraries = $env:R_LIBS
$separator = [IO.Path]::PathSeparator
$env:R_LIBS = $releaseLibrary + $separator + $effectiveLibraries
Rscript --vanilla -e "library(ctgimme); stopifnot(as.character(packageVersion('ctgimme')) == '0.0.12'); stopifnot(identical(sort(getNamespaceExports('ctgimme')), c('ctgimme', 'ctgimme_demo'))); print(ctgimme_demo()[['membership']]); print(citation('ctgimme'))"
$smokeExit = $LASTEXITCODE
$env:R_LIBS = $priorLibraries
if ($smokeExit -ne 0) { throw "Installed archive smoke test failed." }
```

Before CRAN submission, also require clean R-devel and current-release checks
on win-builder and successful Linux, Windows, and macOS GitHub Actions jobs.

## Installed-package ISDSA acceptance

The repository-only scripts under `validation/` exercise the release archive
after it has been installed into an isolated library, without loading
repository source. The simulated fixture and frozen reference CSV files are
included there with checked SHA-256 hashes. The entire directory is excluded
from the CRAN tarball and ordinary `R CMD check` runs.

The exact-support multisubject acceptance compares both ISDSA clone subgroups
in forward and reversed subject order against the frozen 0.0.11 baseline and
generates both fitted RDS files plus eight parameter/transition PNGs. The
exact final archive, installed into a clean library under R 4.6.1/OpenMx
2.22.11 and compared with the R 4.4.2/OpenMx 2.21.13 reference, gave:

- all four fits at status 0;
- ten shared parameters and ten independently initialized likelihood blocks
  per subgroup;
- maximum estimate difference `5.07e-9`;
- maximum standard-error difference `2.96e-6`;
- maximum likelihood difference `4.00e-11`;
- maximum forward/reverse estimate difference `2.83e-10`;
- truth-recovery MAE/RMSE of `0.03584/0.04344` for subgroup 1 and
  `0.01889/0.02413` for subgroup 2.

Run `validation/validate-isdsa-installed.R` after installing the exact final
archive into an isolated library, or trigger the manual GitHub Actions
workflow `Installed ISDSA release acceptance`. The slower
`validation/validate-isdsa-public-installed.R` additionally checks the public
search/orchestration path. The exact final archive completed this gate in
433.28 seconds: it recovered exact membership, group support, and both
subgroup supports; all 20 individual support masks matched the retained
installed 0.0.10 practice run. All 20 individual models converged, 19 matched
truth exactly, and subject 5 retained the same additional finite-sample path
seen in 0.0.8 and 0.0.10 while missing none of the true paths. Both saved
subgroup models contained ten independent subject
likelihoods and one shared ten-parameter vector. Retained evidence is under
`dist/cran-final-0.0.12-h/isdsa-public-installed-validation/`. Do not
reinterpret this single stochastic fixture as proof of frequentist
unbiasedness.

## Git publication

Review and commit the resubmission branch, publish it for CI, and merge only
after its checks succeed:

```sh
git switch cran-resubmission/ctgimme-0.0.12
git add -A
git commit -m "Prepare ctgimme 0.0.12 CRAN resubmission"
git push -u origin cran-resubmission/ctgimme-0.0.12
# Open a pull request and manually wait for all five matrix jobs:
# macOS release, Windows release, and Ubuntu release/devel/oldrel.
git switch main
git merge --ff-only cran-resubmission/ctgimme-0.0.12
git push origin main
# Wait for the main-branch workflow to complete successfully.
git push origin pre-package-0.0.6 v0.0.11
git tag -a v0.0.12 -m "ctgimme 0.0.12"
git push origin v0.0.12
```

Create the GitHub release from `v0.0.12` and attach the checked
`dist/ctgimme_0.0.12.tar.gz`. Do not commit source archives or check
directories to the repository.

Update the GitHub repository About panel from its temporary-repository text:

- Description: `Continuous-time GIMME for group-, subgroup-, and
  individual-level dynamic networks in intensive longitudinal data.`
- Topics: `r`, `r-package`, `time-series`, `continuous-time`,
  `state-space-models`, `gimme`, `intensive-longitudinal-data`.

The repository currently has no branch protection. After this release,
protect `main` and require the five-job R-CMD-check workflow before merging.

Suggested release title: `ctgimme 0.0.12`

Suggested release summary:

> ctgimme 0.0.12 renames the package and primary functions to `ctgimme`, adds
> suppressible progress through `verbose`, replaces the non-executable example
> with a fast self-contained fit, and caps package parallelism at two workers.
> It retains one shared-parameter multisubject likelihood per subgroup: each
> member has its own time series and filter initialization, `time.intervals`
> generates `exp(A * delta)` transition plots, and each subgroup writes one
> `Subgroup_<g>Model.RDS`.

## CRAN resubmission

Upload only the exact checked `ctgimme_0.0.12.tar.gz`. In the submission form,
state that the package was renamed from the rejected `ctsgimme` submission
before its first CRAN publication and summarize each reviewer-requested change.
The detailed response is maintained in `cran-comments.md`.
