## Release

This is ctgimme 0.1.0, the first feature update after the CRAN publication of
0.0.12. The public function signatures and statistical model specification
are unchanged.

## Changes

- The public `cores` argument still defaults to `1`, so examples and ordinary
  default execution remain serial. Explicit user requests are no longer
  subject to a package-wide two-worker ceiling; they are reduced only when
  they exceed the number of subjects. OpenMx uses one thread per R process.
- Warnings raised during subject fits on PSOCK workers are returned to the
  main process and signaled there as ordinary R warnings. Informational
  worker output remains suppressed, and `verbose = FALSE` remains quiet on
  successful runs.
- No fitting, path-selection, multisubject subgroup-model, or saved-output
  behavior was intentionally changed.

Package examples use the serial default. Real process-spawning tests use no
more than two workers; higher requested counts are covered with mocked
worker-pool tests. A separate installed-package smoke test used four workers
outside `R CMD check`.

## Local checks

Windows 11 x64, R 4.4.2:

`R CMD check --as-cran --run-donttest --timings`

0 errors | 0 warnings | 3 notes

1. CRAN incoming reported `Days since last update: 0`, because 0.0.12 was
   published on the day this candidate was prepared.
2. The local system clock could not be verified.
3. Pandoc was unavailable locally, so the top-level README/NEWS Markdown
   check was skipped. Package HTML and PDF manuals were generated and checked
   successfully.

All examples and 213 test expectations passed. The `ctgimme()` example took
0.46 seconds. The checked archive is `ctgimme_0.1.0.tar.gz`, 52,096 bytes,
SHA-256
`C23BC88796F2207FC16F4F7E0CD0A39C00A78D534BE9EDAB1F8AB330921C0005`.

## Installed-package regression validation

The exact archive was installed into an isolated library and compared with
the frozen 0.0.11 multisubject ISDSA reference. All four subgroup fits
converged with ten shared parameters and ten independent subject likelihood
blocks. Maximum differences were `4.44e-16` for estimates, `4.86e-17` for
standard errors, `4.37e-11` for likelihoods, and `5.18e-9` under reversed
subject order. All truth-recovery and PNG/RDS artifact gates passed.

An installed public toy workflow also produced identical membership, group
structure, and individual drift estimates at `cores = 1` and `cores = 4`.
The four-worker run used four distinct PSOCK worker processes.

Current R/R-devel, win-builder, and cross-platform GitHub results will be
added here before CRAN submission. Because 0.0.12 was published today, this
0.1.0 archive should not be submitted immediately unless CRAN specifically
requests a correction.
