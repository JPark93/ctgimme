# ctgimme 0.0.12

- Renamed the package and its complete public identity from `ctsgimme` to
  `ctgimme`. The primary functions are now `ctgimme()` and `ctgimme_demo()`;
  internal helper prefixes, returned-object attributes, help aliases, tests,
  citation metadata, and source archives use the same name. The old function
  names are not exported.
- Added the public `verbose` argument. `verbose = FALSE` suppresses package
  progress messages and captures OpenMx optimizer/progress output, while
  warnings, errors, returned results, and saved artifacts remain available.
  Removed the unconditional console print of subgroup membership.
- Replaced the non-executable wrapped example with a self-contained,
  timed toy fit that writes below the session temporary directory and always
  cleans its output.
- Capped package-managed parallel execution at two workers in all contexts,
  strengthened documented input and control validation, and guaranteed an
  extractable membership object when subgroup detection falls back to one
  group.
- Added a complete, ID-aligned `ctgimme.membership` attribute to every
  successful PAM or legacy clustering return; PAM results also retain their
  direct `membership` element.
- Propagated the master session's library paths to PSOCK workers and exported
  the refactored initial-covariance helpers, preventing two-worker fits from
  failing after an R/library upgrade or during subject-model construction.
- Relayed warning conditions raised by PSOCK subject-fit workers through the
  main R process with `cores = 2`, without changing fitting or selection
  behavior.
- Hardened subject-artifact handling for portable filesystems: identifiers
  must be unique without regard to case, generated prefixes/suffixes are
  decoded literally, wildcard characters are never expanded during cleanup,
  and failed RDS fallback copies now raise an error.
- Retained the 0.0.11 multisubject subgroup-model implementation: one shared
  parameter vector per subgroup, independently initialized subject likelihood
  blocks, requested `exp(A * delta)` plots, and one
  `Subgroup_<g>Model.RDS` artifact.

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
