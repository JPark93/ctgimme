## Test environments

- Windows 11 x64, R 4.4.2 (ucrt)

## R CMD check results

`R CMD check --no-manual` for the ctsgimme 0.0.11 source package on
2026-08-19
completed with `Status: OK` (0 errors, 0 warnings, and 0 notes). The source
archive was created with `R CMD build` under R 4.4.2. A fresh `--as-cran`
check is still required before any CRAN submission.

## Additional validation

The source archive installs into a fresh library. The deterministic package
demo completes successfully. The unit suite contains 158 passing expectations
with no failures, warnings, errors, or skips under the check environment.

Live process-level tests verify that persistent PSOCK workers terminate after
normal completion and error unwinding. They also inject a graceful
`stopCluster()` failure and verify that the package force-terminates the exact
surviving process handles.

The final multisubject implementation was evaluated on the 20-subject,
4,000-row ISDSA fixture using the exact true subgroup supports. All four final
fits (two subgroups in forward and reversed subject order) converged with
status 0, ten shared free parameters, and ten independent subject likelihood
blocks. Estimates were invariant to reversed subject order up to 5.2e-9.
Historical pre-final comparisons against the removed chained estimator are
retained outside the package as development evidence.

The final output validation also regenerates both subgroup parameter graphs,
six requested delta-t transition graphs, and the two subgroup-specific fitted
model RDS files. The graph margins were visually checked to ensure boundary
nodes, loops, titles, and labels are retained.

## Runtime controls

The default is one worker. When `_R_CHECK_LIMIT_CORES_` is active,
package-managed parallelism is limited to two workers, with one OpenMx thread
per worker. Explicit larger requests are allowed outside package checks and
are bounded by the number of subjects; they are not separately capped at the
number of available logical processors. One PSOCK pool is reused across
fitting batches and is stopped on completion or error, with a
process-handle-based forced-shutdown fallback. The unit suite includes one
compact real OpenMx subgroup fit that verifies the parameter plot, two
requested delta-t transition plots, and subgroup-specific fitted-model RDS.
