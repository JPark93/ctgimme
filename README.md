# ctgimme

[![R-CMD-check](https://github.com/JPark93/ctgimme/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JPark93/ctgimme/actions/workflows/R-CMD-check.yaml)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)

`ctgimme` estimates group-, subgroup-, and individual-level continuous-time
dynamic networks from intensive longitudinal data using continuous-time
subgrouping GIMME (C-TSGIMME). It combines continuous-time state-space models
with iterative searches for temporal paths supported across the sample,
shared within recovered membership groups, or specific to individual
subjects.

The method is described in:

> Park, J. J., Fisher, Z. F., Hunter, M. D., Shenk, C., Russell, M.,
> Molenaar, P. C. M., and Chow, S.-M. (2025). Unsupervised model construction
> in continuous-time. *Structural Equation Modeling: A Multidisciplinary
> Journal, 32*(3), 377--399.
> <https://doi.org/10.1080/10705511.2024.2429544>

## Installation

Install the current release directly from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("JPark93/ctgimme", dependencies = TRUE)
```

Install the 0.1.0 source archive while also resolving dependencies:

```r
remotes::install_local(
  "path/to/ctgimme_0.1.0.tar.gz",
  dependencies = TRUE
)
```

If the dependencies are already installed, base R can install the archive:

```r
install.packages(
  "path/to/ctgimme_0.1.0.tar.gz",
  repos = NULL,
  type = "source"
)
```

The imported dependencies, including OpenMx, are listed in `DESCRIPTION`.
Package `expm` is required when `subgroup.model = TRUE`; `igraph` and `nloptr`
support the optional legacy subgrouping workflow. Building or checking from a
source checkout on Windows may require the Rtools version appropriate for the
installed version of R.

## Data requirements

The primary function, `ctgimme()`, expects a data frame in long format with
one row per observation. The `id` and `time` arguments name the subject and
observation-time columns, and `varnames` names the modeled variables. These
columns may appear in any position. Identifier and time values must be
nonmissing; observation times and nonmissing modeled values must be finite and
numeric. Rows are not sorted internally, so they must already be ordered by
time within subject. Missing modeled values are passed to OpenMx's raw-data
likelihood.

The order of `varnames` determines the order of variables in the fitted drift
and noise matrices. Subject identifiers are also used in intermediate
filenames and therefore must be valid portable filename components. The
identifier values must also remain unique when compared without case, because
common Windows and macOS filesystems are case-insensitive. The selected time
column is copied internally to a reserved column named `Time`; do not use
`Time` as the name of a modeled process variable.

## Basic use

```r
library(ctgimme)

result <- ctgimme(
  varnames = c("x1", "x2", "x3"),
  dataframe = my_long_data,
  id = "subject_id",
  time = "observation_time",
  subgroup.method = "pam",
  sub.sig.thrsh = 0.55,
  max.subgroups = 4,
  cores = 1,
  directory = "ctgimme-output",
  verbose = TRUE
)
```

Set `verbose = FALSE` to suppress package progress messages and OpenMx
optimizer reporting. Warnings and errors remain available through standard R
condition handling, and the returned object and written artifacts are the same
under either setting.

## Search and subgroup controls

`sig.thrsh = 0.55` requires a candidate group path to have significant
modification-index support from at least 55% of subjects with usable results.
`Galpha = 0.05` supplies the base group-level significance threshold. The
analogous subgroup defaults are `sub.sig.thrsh = 1` and `S.Galpha = 0.05`, and
the individual-level base threshold is `Ialpha = 0.01`. With the default
`ben.hoch = TRUE`, these alpha values seed an adapted, progressively stricter
threshold sequence for path additions; `FALSE` uses the corresponding alpha at
every step. A group- or membership-group shared search adds no path when usable
modification-index results are available for no more than half of the subjects
requested for that search.

The default `sub.sig.thrsh = 1` bypasses data-driven subgroup detection and
assigns all subjects to membership group 1. It does not skip the subsequent
within-group shared-path search or individual-model fits. Set it below one to
activate `subgroup.method`. The default method, `"pam"`, screens recurrent
signed evidence for eligible non-group paths, computes mean Manhattan
distances, and selects the candidate with the largest average silhouette
width. Active PAM subgrouping requires `max.subgroups` to be an integer of at
least two; the actual upper candidate count is also limited to one less than
the number of subjects with usable features. Its recurrence screen uses fixed
`0.05` subject-evidence and BH-adjusted recurrence thresholds, independently
of `ben.hoch`.

The alternative `subgroup.method = "legacy"` uses a weighted similarity graph
and Walktrap community detection and requires `igraph`. With `conduct = TRUE`,
the default, its graph cutoff is optimized with `nloptr`; set `conduct = FALSE`
to use the unoptimized similarity weights shifted to start at zero.

Set `scale.data = TRUE` to standardize each modeled variable separately within
subject before all estimation stages. The default is `FALSE`; constant
within-subject columns are centered to zero when scaling is requested.

## Noise variances

By default, measurement-error and process-noise variances are fixed at the
values supplied through `ME.var` and `PE.var`: `ME.var` defaults to `1e-8`,
whereas `PE.var = NULL` produces an identity matrix. Set `ME.free = TRUE` or
`PE.free = TRUE` to estimate every diagonal variance, using the corresponding
variance argument as its starting value. Logical vectors allow selective
estimation in the order given by `varnames`:

```r
result <- ctgimme(
  varnames = c("x1", "x2", "x3"),
  dataframe = my_long_data,
  id = "subject_id",
  time = "observation_time",
  directory = "ctgimme-output",
  ME.var = 0.05,
  ME.free = FALSE,
  PE.var = 1,
  PE.free = c(TRUE, TRUE, FALSE)
)
```

Only diagonal measurement-error and process-noise variances are supported;
off-diagonal covariances remain fixed to zero. Estimating both sets of
variances can require substantially more within-person information than
estimating process noise while fixing measurement error. Estimated
process-noise and measurement-error variances are stored in the `Q` and `R`
matrices, respectively, of each saved OpenMx model.

## Optional joint membership-group models

The default `subgroup.model = FALSE` skips only the additional joint,
parameterized model for each membership group; subgroup-level structural
searches and individual fits still run. Set it to `TRUE` to fit one
shared-parameter model per subgroup and write its coefficient plot,
discrete-time transition plots, and fitted-model RDS file. This workflow
requires `expm`.

```r
result <- ctgimme(
  ...,
  subgroup.model = TRUE,
  time.intervals = c(0.25, 1, 2)
)
```

The subgroup model is estimated from the summed likelihoods of its members.
Each member has an independently initialized continuous-time
state-space filter, so no state is propagated across subject boundaries and no
trajectories are concatenated. The saved result is one fitted OpenMx model with
one top-level `A`, `Q`, and `R`; its subject children are likelihood blocks,
not separately estimated parameter models. This creates one fitted joint model
for each subgroup, not one fitted parameter model per member. Each subject's
observation times are rebased to its first observation, preserving all
within-person elapsed intervals while removing arbitrary calendar offsets.

`time.intervals` remains the only subgroup-model time control. For every
requested nonnegative delta t, ctgimme computes `exp(A * delta_t)` and writes
the discrete-time transition plot. For subgroup `g`, the files are written to
`Models/Subgroup <g>/` with these names:

- `Subgroup <g> Params.png`;
- `Subgroup <g> Delta_t = <delta>.png` for every requested delta t; and
- `Subgroup_<g>Model.RDS`, containing the single fitted joint model.

## OpenMx and RStudio compatibility

All internal `mxTryHard()` calls pass through one package wrapper. With
`verbose = TRUE`, package stage messages and ordinary optimizer reporting from
fits in the main R process are shown; PSOCK workers do not forward raw console
output. With `verbose = FALSE`, package messages and main-process OpenMx output
are suppressed so a quiet run remains quiet. Subgroup parameter plots also
supply edge labels as a square matrix, as required by `qgraph` for the fitted
drift matrix. Parameter and transition plots use expanded margins so boundary
nodes, loops, titles, and labels are not clipped in the saved PNG files.

For a first run on a new computer, use `cores = 1` and a new, empty output
directory. After that run succeeds, `cores > 1` can reduce elapsed time, at the
cost of the memory used by each additional R/OpenMx worker. Explicit requests
are reduced only when they exceed the number of subjects; the package imposes
no separate worker maximum. Choose a value supported by the machine or
computing allocation being used.

## Results and output files

When enabled subgroup detection succeeds, the function returns the selected
PAM or Walktrap clustering object. Both forms carry the complete, ID-aligned
vector in `attr(result, "ctgimme.membership")`; PAM results additionally expose
the same vector as `result$membership`. If PAM has fewer than three usable subject
model/MI pairs or no recurrent eligible features, or legacy detection has fewer
than two usable pairs, detection falls back to one membership group and the
function invisibly returns a result list containing the complete membership,
group structure, output directory, and subgroup diagnostics. A successfully
computed legacy solution may itself contain one community and still returns
its communities object.

When `sub.sig.thrsh = 1`, the same result-list form includes the group
structure, complete membership vector, output directory, subgroup-detection
record, and completion message. Thus every successful return directly exposes
the complete membership. Every successfully completed run
also saves it to `subgroup_membership.csv` and `subgroup_detection.rds`.
Subjects without usable model/MI artifacts are assigned to group 1 in these
complete outputs.

The output directory includes:

- `GStruc.RDS`, containing the group-level structure;
- `SGStruc*.RDS`, containing the shared structure for each membership group,
  including the single group used when detection is disabled or falls back;
- `subgroup_membership.csv`, mapping every subject identifier to a subgroup;
- `Subgroups Plot.png`, showing a subject-distance map when diagnostic
  distances are available and a complete membership roster otherwise;
- `subgroup_detection.rds`, containing subgroup diagnostics; and
- successfully fitted final subject models and plots in `Models/Individuals`.

Successful PAM detection also writes `pam_silhouette_selection.png`; successful
legacy detection writes `walktrap_community_plot.png`. With
`subgroup.model = TRUE`, each successfully fitted membership-group model is
written as `Subgroup_<g>Model.RDS`, together with `Subgroup <g> Params.png` and
one `Subgroup <g> Delta_t = <delta>.png` for every requested delta t. By
default, `keep.intermediate = FALSE` removes temporary per-subject
`Model_<id>.RDS` and
`MI_<id>.RDS` search files after successful completion; final individual and
pooled membership-group models are retained.

For PAM subgrouping, the plot uses classical multidimensional scaling of the
same mean Manhattan distances used for clustering. Mapped subjects that are
closer in the plot have more similar recurrent-evidence profiles. Colors and
enclosing hulls show subgroup membership, while nearest-neighbor edges aid
visual interpretation. The two-dimensional display is an approximation; its
reported goodness-of-fit summarizes how well it represents the full distance
matrix. Subjects lacking a diagnostic distance are named on the plot but are
not placed on the map.

## Parallel execution

The `cores` argument defaults to one. Explicit requests are bounded by the
number of subjects only; the package imposes no additional worker maximum and
does not override the user's resource choice. OpenMx is set to one thread in
each R process during the run, and the calling process's previous OpenMx
thread setting is restored on exit. Package examples use the serial default,
and package tests never start more than two workers.

For `cores > 1`, `ctgimme` creates one PSOCK worker pool and reuses it for the
group, subgroup, and individual fitting batches. Workers load OpenMx and
qgraph once, wait idle while the main R process pools modification indices and
selects the next path or subgroup, and then receive the next fitting batch.
Warnings raised by worker fits are relayed through the main R process and
remain suppressible with standard R condition handling. Informational worker
output is not forwarded. The pool is stopped when the analysis finishes or
exits with an error.

## Quick demonstration

`ctgimme_demo()` provides a deterministic, lightweight demonstration of the
feature-selection and subgrouping logic without fitting OpenMx models:

```r
demo <- ctgimme_demo()
demo$membership
demo$candidates
```

The demonstration avoids OpenMx fitting and deterministically recovers its
known two-subgroup partition.

## Citation

Use the package citation and its methodological reference with:

```r
citation("ctgimme")
```

## Development

See the
[contributing guide](https://github.com/JPark93/ctgimme/blob/main/CONTRIBUTING.md)
for the test and pull-request workflow.
Release notes are maintained in [NEWS.md](NEWS.md).

## Authors and license

`ctgimme` is authored by Jonathan J. Park
(<imJPark@UCDavis.edu>) and Nathan Xin Mills
(<nxmills@berkeley.edu>). The package is distributed under the Apache License
2.0.
