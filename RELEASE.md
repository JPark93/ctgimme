# Publishing ctgimme 0.1.0

Version 0.1.0 is the first feature release after the CRAN publication of
0.0.12. It preserves the public API and statistical fitting behavior while
making two parallel-execution changes:

- explicit `cores` requests are no longer subject to a package-wide
  two-worker ceiling (the default remains `cores = 1`); and
- warnings raised during subject fits on PSOCK workers are relayed through
  the main R process as ordinary, suppressible R warnings.

OpenMx still uses one thread in each R process. Requested worker counts are
reduced only when they exceed the number of subjects.

## Checked archive

The checked source archive is `dist/ctgimme_0.1.0.tar.gz` (52,096 bytes),
SHA-256
`C23BC88796F2207FC16F4F7E0CD0A39C00A78D534BE9EDAB1F8AB330921C0005`.

On Windows 11 with R 4.4.2,
`R CMD check --as-cran --run-donttest --timings` completed with 0 errors,
0 warnings, and three local notes: the same-day CRAN update status, inability
to verify the local system clock, and unavailable Pandoc for checking the
top-level Markdown files. Package code, installation, examples, tests, PDF
and HTML manuals, and documentation checks all passed. The suite reported
213 passing expectations; the primary example elapsed time was 0.46 seconds.

Retained evidence is under `dist/ctgimme-0.1.0-final/`. Do not rebuild or
overwrite the accepted archive after publishing its hash.

## Installed-package acceptance

The exact archive was installed into a fresh library and validated against
the retained ISDSA clone fixture and frozen 0.0.11 multisubject estimates.
All convergence, architecture, parameter/standard-error/likelihood parity,
subject-order invariance, truth-recovery, demo, and PNG/RDS artifact gates
passed. The run took 63.92 seconds and reported:

- maximum estimate difference from the reference: `4.44e-16`;
- maximum standard-error difference: `4.86e-17`;
- maximum likelihood difference: `4.37e-11`; and
- maximum forward/reverse estimate difference: `5.18e-9`.

An additional installed-package smoke test created four distinct PSOCK
workers and compared the public toy workflow at `cores = 1` and `cores = 4`.
Membership and group structure were identical, and the maximum difference
among fitted individual drift estimates was zero.

The slower full public ISDSA result for 0.0.12 remains recorded in
`validation/ACCEPTANCE-0.0.12.md`. It is historical evidence, not a claim
that the full public gate was rerun for 0.1.0.

## Reproducing local checks

Run from a clean release commit with current R, LaTeX, and Pandoc available:

```powershell
Rscript --vanilla -e "devtools::test('.', stop_on_failure = TRUE)"
Rscript --vanilla -e "spelling::spell_check_package('.')"
Rscript --vanilla -e "urlchecker::url_check('.')"

New-Item -ItemType Directory -Force dist/ctgimme-0.1.0-check | Out-Null
Push-Location dist/ctgimme-0.1.0-check
R.exe CMD build ../.. --no-manual
R.exe CMD check --as-cran --run-donttest --timings ctgimme_0.1.0.tar.gz
Pop-Location
```

Install and smoke-test the accepted archive in a new library:

```powershell
$archive = (Resolve-Path "dist/ctgimme_0.1.0.tar.gz").Path
$expectedSha = "C23BC88796F2207FC16F4F7E0CD0A39C00A78D534BE9EDAB1F8AB330921C0005"
if ((Get-FileHash -Algorithm SHA256 -LiteralPath $archive).Hash -ne $expectedSha) {
  throw "Release archive SHA-256 mismatch."
}
$library = Join-Path $env:TEMP ("ctgimme-0.1.0-" + [guid]::NewGuid())
New-Item -ItemType Directory -Path $library | Out-Null
R.exe CMD INSTALL "--library=$library" $archive
if ($LASTEXITCODE -ne 0) { throw "R CMD INSTALL failed." }
$effectiveLibraries = Rscript --vanilla -e 'cat(paste(.libPaths(), collapse = .Platform$path.sep))'
$priorLibraries = $env:R_LIBS
$env:R_LIBS = $library + [IO.Path]::PathSeparator + $effectiveLibraries
Rscript --vanilla -e "library(ctgimme); stopifnot(as.character(packageVersion('ctgimme')) == '0.1.0'); stopifnot(identical(sort(getNamespaceExports('ctgimme')), c('ctgimme', 'ctgimme_demo'))); print(ctgimme_demo()[['membership']])"
$smokeExit = $LASTEXITCODE
$env:R_LIBS = $priorLibraries
if ($smokeExit -ne 0) { throw "Installed archive smoke test failed." }
```

The full ISDSA command is documented in `validation/README.md`.

## Git publication

Push the release branch and wait for all five R-CMD-check matrix jobs before
merging. The repository currently has no required-check branch protection,
so this is a manual gate.

```sh
git switch release/0.1.0
git push -u origin release/0.1.0
# Open a pull request and wait for macOS release, Windows release,
# and Ubuntu release/devel/oldrel.
git switch main
git merge --ff-only release/0.1.0
git push origin main
# Wait for the main-branch checks, then tag the exact accepted commit.
git tag -a v0.1.0 -m "ctgimme 0.1.0"
git push origin v0.1.0
```

Create the GitHub release from `v0.1.0` and attach only
`dist/ctgimme_0.1.0.tar.gz`. Do not commit source archives, installed
libraries, or check directories.

Suggested release summary:

> ctgimme 0.1.0 keeps serial execution as the default while honoring larger
> explicit worker requests up to the number of subjects. OpenMx remains at
> one thread per R process. Parallel subject-fit warnings are now relayed to
> the caller without changing fitting, selection, or saved-output behavior.

## CRAN timing

CRAN published 0.0.12 on the same day this candidate was prepared. Do not
submit 0.1.0 immediately unless CRAN requests a correction. First allow the
0.0.12 check page to settle for at least 48 hours, then complete current R,
R-devel/win-builder, and cross-platform GitHub checks on the exact 0.1.0
commit. Update `cran-comments.md` with those final external results before
submission.
