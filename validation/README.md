# Installed-package ISDSA acceptance test

`validate-isdsa-installed.R` verifies the installed source archive rather than
loading the repository source. It reproduces the truth-controlled ISDSA clone
comparison used to validate the multisubject subgroup likelihood in 0.0.11.
`validate-isdsa-public-installed.R` is the slower orchestration gate: it runs
the installed public `ctgimme()` workflow and checks membership, group,
subgroup, individual, joint subgroup-model, quiet-output, cleanup, and artifact
contracts against the fixture truth.

The checked 0.0.12 outcomes are summarized in
[`ACCEPTANCE-0.0.12.md`](ACCEPTANCE-0.0.12.md).

The repository includes the simulated fixture and frozen 0.0.11 numerical
reference under `validation/fixtures/` and `validation/reference-0.0.11/`.
Their SHA-256 hashes are checked before estimation. You may override those
defaults with these arguments:

- `--fixture`: `isdsa_jpark26_fixture.rds` (expected SHA-256
  `5ADD83906470012F9DCE2B0F89098683B02370A42E698EDD018C9A8B24412985`);
- `--reference`: the prior `analysis/demo-multisubject` directory containing
  `fit_summary.csv`, `parameter_estimates.csv`, and
  `subject_order_invariance.csv`;
- `--output`: a required new directory for 0.1.0 results and generated
  artifacts.

Install the exact tarball into an isolated library and put that library first
on the child R process's library path before running the harness. The
acceptance fit also requires the suggested package `expm`; install it once if
it is not already available. From the repository root on Windows PowerShell:

```powershell
$archive = (Resolve-Path "dist/ctgimme_0.1.0.tar.gz").Path
$expectedArchiveSha256 = "C23BC88796F2207FC16F4F7E0CD0A39C00A78D534BE9EDAB1F8AB330921C0005"
if ((Get-FileHash -Algorithm SHA256 -LiteralPath $archive).Hash -ne $expectedArchiveSha256) {
  throw "Release archive SHA-256 mismatch."
}
$runId = [guid]::NewGuid().ToString("N")
$library = Join-Path $env:TEMP ("ctgimme-0.1.0-library-" + $runId)
$output = Join-Path $env:TEMP ("ctgimme-0.1.0-output-" + $runId)
New-Item -ItemType Directory -Path $library -Force | Out-Null
Rscript --vanilla -e "if (!requireNamespace('expm', quietly = TRUE)) install.packages('expm')"
if ($LASTEXITCODE -ne 0) { throw "Installing acceptance dependencies failed." }
R.exe CMD INSTALL "--library=$library" $archive
if ($LASTEXITCODE -ne 0) { throw "R CMD INSTALL failed." }

$effectiveLibraries = Rscript --vanilla -e 'cat(paste(.libPaths(), collapse = .Platform$path.sep))'
if ($LASTEXITCODE -ne 0) { throw "Reading the R library path failed." }
$priorLibraries = $env:R_LIBS
$separator = [IO.Path]::PathSeparator
$env:R_LIBS = $library + $separator + $effectiveLibraries
Rscript --vanilla validation/validate-isdsa-installed.R `
  --output=$output
$acceptanceExit = $LASTEXITCODE
$env:R_LIBS = $priorLibraries
if ($acceptanceExit -ne 0) { throw "Installed ISDSA acceptance failed." }
```

The script records and reports the exact installed package path. Use a new
library and output directory for every release-candidate run.

Acceptance requires convergence, the expected shared-parameter/multisubject
architecture, numerical parity with the prior validated estimates, subject-
order invariance, truth-recovery bounds, the built-in demo partition, and the
complete subgroup PNG/RDS artifact contract. The public gate additionally
requires all 20 individual support masks to match the retained installed
0.0.10 practice run. PNG and RDS files are checked
semantically rather than byte-for-byte because graphics and serialized model
metadata vary across R, OpenMx, and platform versions.

The full public gate takes several minutes and is intentionally excluded from
ordinary `R CMD check` and pull-request CI. Run it manually for release
candidates with the same `--fixture` and `--output` arguments (it does not need
`--reference`).
