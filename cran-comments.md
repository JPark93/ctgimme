## Resubmission

This is a resubmission of the previously rejected `ctsgimme` package. The
package was not published on CRAN under that name. Before resubmission, the
package, namespace, primary function, demo function, documentation aliases,
returned-object attributes, tests, citation metadata, and source archive were
renamed consistently to `ctgimme`. The resubmission version is 0.0.12.

## Response to reviewer comments

### Examples

- Removed the `\dontrun{}` wrapper.
- Replaced the incomplete example, which referred to an undefined object,
  with a self-contained and deterministic toy fit.
- The example is executable, is therefore left unwrapped, writes only beneath
  the R session temporary directory, and cleans its files on exit.
- The exact submitted archive completed the `ctgimme()` example in 0.68
  seconds during `R CMD check --as-cran --run-donttest`.
- The package contains no `\dontrun{}` or `\donttest{}` examples.

### Console output

- Removed the unconditional `print()` call that displayed subgroup
  membership. Membership is returned by `ctgimme()` and remains available in
  the saved membership artifact.
- Added a scalar logical `verbose` argument at the end of the public function
  signature. With `verbose = FALSE`, package progress messages and OpenMx
  optimizer/progress output are suppressed. Results and saved artifacts are
  unchanged.
- User-facing progress is emitted with `message()` only when `verbose = TRUE`.
  Warnings and errors remain proper, suppressible R conditions.
- Tests exercise both verbose modes, including real serial and two-worker
  fits, and verify that `verbose = FALSE` produces no standard output or
  messages.

## Additional changes and checks

- Package-managed parallelism is capped at two workers in every environment,
  with one OpenMx thread per worker.
- Added validation for identifiers, time and variable columns, controls, and
  output paths; hardened temporary-file cleanup and portable subject artifact
  names.
- Added worker-library propagation and all required PSOCK exports, with a real
  two-worker regression test.
- Retained the multisubject subgroup estimator introduced in the preceding
  0.0.11 package version: each subgroup has one shared parameter vector and one
  independently initialized likelihood block per subject. Requested
  `exp(A * delta)` plots and one fitted `Subgroup_<g>Model.RDS` are produced
  per subgroup.

The exact release archive also passed an installed-package acceptance test on
the 20-subject, 4,000-row ISDSA clone fixture. All four subgroup fits (two
subgroups in forward and reversed subject order) converged with ten shared
parameters and ten independent subject likelihood blocks. Relative to the
frozen 0.0.11 multisubject reference, the maximum estimate difference was
5.07e-9, the maximum standard-error difference was 2.96e-6, the maximum
likelihood difference was 4.00e-11, and the maximum forward/reverse estimate
difference was 2.83e-10. This is a finite-fixture
regression test and is not presented as evidence of frequentist unbiasedness.

The installed public `ctgimme()` workflow was also run on the complete fixture
with `verbose = FALSE`. It produced no console output or messages, recovered
the exact membership plus group and subgroup supports, converged for all 20
individual fits, and wrote two joint ten-subject subgroup models with all
requested RDS and transition-plot artifacts.

## Test environment

- Windows 11 x64 (build 26200)
- R 4.6.1 (2026-06-24 ucrt)
- `R CMD check --as-cran --run-donttest --timings`

## R CMD check results

0 errors | 0 warnings | 2 notes

1. `New submission` -- expected because `ctgimme` has not previously been
   published on CRAN.
2. Local HTML validation was skipped because an HTML Tidy executable was not
   installed. The HTML manual was generated successfully; this note is local
   environment configuration rather than a package issue.

The PDF manual, examples, demos, and tests completed successfully. The test
suite reported 190 passing expectations with no failures, warnings, errors,
or skips.

The checked archive is `ctgimme_0.0.12.tar.gz`, SHA-256
`485E698808387C8FB15E30EC9B464878BCF88B82BEF710ED2D756CE955AF1ABA`.
