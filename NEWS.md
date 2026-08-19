# ctsgimme 0.0.11

- Converted the historical single-file repository into an installable R
  package with generated help, tests, citation metadata, and standard GitHub
  installation. The 0.0.6 monolithic script and generated example outputs
  remain recoverable from the `pre-package-0.0.6` Git tag.
- Replaced the previous concatenated optional subgroup fit with one
  shared-parameter multisubject likelihood. Each subject receives an
  independent continuous-time state-space filter while the subgroup has one
  shared `A`/`Q`/`R` specification whose free entries are jointly optimized.
  Subject times are rebased locally to preserve elapsed intervals without
  making the likelihood depend on arbitrary start offsets.
- Removed the chained estimator and its public selector and concatenation
  controls: `subgroup.model.method`, `subgroup.time.mode`,
  `measurement.schedule`, `cycle.interval`, `insert.na.rows`, and
  `subject.gap`. This intentionally changes calls that supplied those 0.0.10
  arguments, including positional calls after `time.intervals`.
- Retained `time.intervals` for producing `exp(A * delta)` discrete-time
  transition plots. Every successfully fitted subgroup also writes its
  parameter plot and `Subgroup_<g>Model.RDS` joint fitted model.
- Expanded the margins on subgroup parameter and transition graphs so outer
  nodes, self-loops, titles, and edge labels remain visible in saved PNGs.
- Added an optional internal `P0.values` override and factored the existing
  empirical initial-covariance calculation into reusable helpers.

# ctsgimme 0.0.10

- Reused one initialized PSOCK worker pool across group, subgroup, and
  individual fitting batches, avoiding repeated R-process and package startup
  costs during a single analysis.
- Ensured package-owned worker pools are stopped after successful completion
  and during error unwinding.

# ctsgimme 0.0.9

- Added `ME.free` and `PE.free` options for estimating selected diagonal
  measurement-error and process-noise variances. Both default to `FALSE`, and
  off-diagonal noise covariances remain fixed to zero.
- Prevented OpenMx's interactive `imxReportProgress` callback from being used
  by internal `mxTryHard()` fits. This applies to group, subgroup, and
  individual estimation paths.
- Fixed subgroup parameter plotting by passing the drift edge labels to
  `qgraph` as a square matrix.
- Expanded documentation and tests for the diagonal-noise options and the
  RStudio/OpenMx compatibility fixes.
