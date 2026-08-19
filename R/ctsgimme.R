#' Estimate Continuous-Time Group, Subgroup, and Individual Networks
#'
#' Fits continuous-time state-space models for multiple individuals and uses
#' iterative modification-index searches to identify paths shared at the
#' group, subgroup, and individual levels. When data-driven subgroup detection
#' is enabled, the default method uses recurrent-evidence feature screening
#' followed by partitioning around medoids (PAM).
#'
#' @param varnames Nonempty character vector naming distinct process-variable
#'   columns in `dataframe`. These columns must be numeric. Their order defines
#'   the row and column order of the drift and noise matrices and the order of
#'   vector-valued `ME.free` and `PE.free` specifications. This argument is
#'   required.
#' @param dataframe Data frame in long format, with one row per measurement
#'   occasion and columns named by `id`, `time`, and `varnames`. Identifier and
#'   time values must be nonmissing, time values and observed process values
#'   must be finite and numeric, and rows must already be ordered by time
#'   within subject. Missing process measurements are passed to OpenMx's
#'   raw-data likelihood. This argument is required.
#' @param id Character scalar naming the subject-identifier column in
#'   `dataframe`. Each distinct value defines one subject. Values are converted
#'   to character for matching, membership labels, and intermediate filenames,
#'   so use nonmissing identifiers that are safe in filenames. This argument
#'   is required.
#' @param time Character scalar naming the numeric observation-time column in
#'   `dataframe`. Its units determine the scale of continuous-time parameters
#'   fit to the observed-time data. The multisubject subgroup model preserves
#'   the supplied elapsed intervals independently within each subject and
#'   rebases each person's first observation to local time zero. These values
#'   also define the units of `time.intervals`. The selected column is copied
#'   internally to a reserved column named `Time`, so `Time` should not name a
#'   process variable. This argument is required.
#' @param cores Positive integer giving the requested number of parallel PSOCK
#'   workers. The default is `1`, which runs without a worker pool. Requests
#'   are reduced to the number of subjects, but are not otherwise capped during
#'   ordinary use. When `_R_CHECK_LIMIT_CORES_` is active during a package
#'   check, they are additionally capped at two. One pool is reused across
#'   fitting stages, and OpenMx uses one thread in each R process.
#' @param directory Character scalar giving the output directory for model,
#'   modification-index, plot, membership, and diagnostic artifacts. It is
#'   created recursively when necessary. Existing files with package-defined
#'   names may be replaced or cleaned up, so a new or dedicated directory is
#'   recommended. This argument is required.
#' @param sig.thrsh Numeric scalar in `[0, 1]` giving the minimum support
#'   proportion for group-level paths. A candidate is added when at least this
#'   proportion of subjects with usable modification-index results support it
#'   at the current `Galpha` threshold; the proportion also governs
#'   low-support pruning. No shared path is added unless more than half of the
#'   requested subjects have usable modification-index results. The default is
#'   `0.55`.
#' @param sub.sig.thrsh Numeric scalar in `[0, 1]` with two roles. Values below
#'   `1` enable the detector selected by `subgroup.method` and set the support
#'   proportion for paths added and retained within each recovered membership
#'   group. The default, `1`, bypasses data-driven detection and assigns every
#'   subject to membership group 1; it does not skip the subsequent
#'   within-group shared-path search or individual-model fitting.
#' @param Galpha Numeric scalar between zero and one giving the base
#'   significance level for modification-index evidence during group-level
#'   path additions. The default is `0.05`; `ben.hoch` determines whether it is
#'   used directly or as the start of an adapted threshold sequence.
#' @param ben.hoch Logical scalar. If `TRUE` (the default), use the adapted
#'   Benjamini-Hochberg threshold sequence during group-, subgroup-, and
#'   individual-level path additions. If `FALSE`, use the corresponding
#'   `Galpha`, `S.Galpha`, or `Ialpha` value at every addition step. This option
#'   does not control the separate, fixed BH correction used to screen
#'   recurrent features for PAM subgrouping.
#' @param S.Galpha Numeric scalar between zero and one giving the base
#'   significance level for modification-index evidence during subgroup-level
#'   path additions. The default is `0.05`; its stepwise use is controlled by
#'   `ben.hoch`.
#' @param Ialpha Numeric scalar between zero and one giving the base
#'   significance level for individual-level path additions and the direct
#'   significance cutoff for individual-level pruning. The default is `0.01`;
#'   `ben.hoch` controls its stepwise use during additions.
#' @param ME.var Finite, nonnegative measurement-error variances, supplied as a
#'   numeric scalar or a diagonal numeric matrix with one row and column per
#'   process variable. A scalar is repeated down the diagonal. The default is
#'   `1e-8`. Entries selected by `ME.free` are starting values and must be at
#'   least `1e-8`; other entries remain fixed. These values populate the
#'   OpenMx `R` matrix.
#' @param PE.var Finite, nonnegative process-noise (innovation) variances,
#'   supplied as `NULL`, a numeric scalar, or a diagonal numeric matrix with one
#'   row and column per process variable. A scalar is repeated down the
#'   diagonal; the default, `NULL`, uses an identity matrix. Entries selected
#'   by `PE.free` are starting values and must be at least `1e-5`; other entries
#'   remain fixed. These values populate the OpenMx `Q` matrix.
#' @param ME.free Logical scalar, vector with one entry per process variable,
#'   or square logical matrix with one row and column per process variable,
#'   selecting the diagonal measurement-error variances to estimate. A scalar
#'   is repeated down the diagonal; a matrix must have `FALSE` off the
#'   diagonal. The default is `FALSE`, which fixes all measurement-error
#'   variances at `ME.var`.
#' @param PE.free Logical scalar, vector with one entry per process variable,
#'   or square logical matrix with one row and column per process variable,
#'   selecting the diagonal process-noise variances to estimate. A scalar is
#'   repeated down the diagonal; a matrix must have `FALSE` off the diagonal.
#'   The default is `FALSE`, which fixes all process-noise variances at
#'   `PE.var`.
#' @param subgroup.model Logical scalar. If `TRUE`, fit one additional joint
#'   continuous-time model for each recovered membership group using that
#'   group's final drift structure. Every member contributes an independent
#'   state-space likelihood and fresh filter initialization. The subgroup has
#'   one shared `A`/`Q`/`R` specification, whose free entries are optimized
#'   jointly once. The fit writes a
#'   parameter plot, the requested discrete-time transition plots, and one
#'   subgroup-specific fitted-model RDS file. No subject series are chained or
#'   separated by synthetic rows. The default is `FALSE`; subgroup structural
#'   searches and individual-model fits are still performed. This option
#'   requires the suggested package `expm`.
#' @param time.intervals Nonempty numeric vector of nonnegative elapsed-time
#'   values in the same units as `time`. When `subgroup.model = TRUE`, every
#'   value `delta` produces a plot of the fitted discrete-time transition
#'   matrix `exp(A * delta)` for each membership group. The default is `1`;
#'   this argument is otherwise ignored.
#' @param keep.intermediate Logical scalar controlling end-of-run cleanup. If
#'   `FALSE` (the default), temporary per-subject `Model_<id>.RDS` and
#'   `MI_<id>.RDS` files from the group- and subgroup-level shared-path
#'   searches are removed after a successful run. If `TRUE`, the final set of
#'   those files is retained; iterations may still replace earlier versions.
#'   Final individual models, fitted subgroup models, and summary artifacts
#'   are retained under either setting.
#' @param conduct Logical scalar used only when `subgroup.method = "legacy"`
#'   and subgroup detection is enabled. If `TRUE` (the default), `nloptr`
#'   selects a quantile cutoff for the subject-similarity graph using the
#'   legacy conductance objective before Walktrap clustering. If `FALSE`,
#'   Walktrap uses the unoptimized similarity weights after shifting their
#'   minimum to zero. `TRUE` requires the suggested package `nloptr`. This
#'   argument is ignored by the default `"pam"` method and when
#'   `sub.sig.thrsh = 1`.
#' @param subgroup.method Character scalar selecting the data-driven subgroup
#'   detector. The default, `"pam"`, screens signed, non-group path evidence for
#'   recurrence, computes mean Manhattan distances, and chooses a
#'   partitioning-around-medoids solution by average silhouette width.
#'   `"legacy"` constructs a weighted subject-similarity graph and applies
#'   Walktrap community detection; it requires the suggested package `igraph`.
#'   PAM's recurrence screen uses a fixed `0.05` subject-level evidence
#'   threshold and retains features with a BH-adjusted recurrence value no
#'   greater than `0.05`. The choice is inactive when `sub.sig.thrsh = 1`.
#' @param max.subgroups `NULL` or a single integer. When PAM subgrouping is
#'   active, an integer at least `2` is required and gives the requested maximum
#'   number of subgroups. Candidate counts run from two through the smaller of
#'   `max.subgroups` and one less than the number of subjects with usable
#'   subgroup features. This argument is ignored for `"legacy"` subgrouping
#'   and when `sub.sig.thrsh = 1`.
#' @param scale.data Logical scalar. If `TRUE`, center and standardize every
#'   process variable separately within each subject before all estimation
#'   stages; within-subject constant columns are centered to zero. The default
#'   is `FALSE`, which preserves the supplied scale.
#'
#' @return If `sub.sig.thrsh = 1`, an invisible list with elements `message`,
#'   `G.DRIFT`, `membership`, and `directory`. With enabled subgroup detection,
#'   a successful PAM run returns a `cluster::pam` object augmented with a
#'   complete, ID-aligned `membership` element; a successful legacy run returns
#'   an `igraph` communities object. Both successful clustering objects carry
#'   a `ctsgimme.membership.artifacts` attribute. If enabled detection falls
#'   back to a single group because too few subject model/modification-index
#'   pairs are usable or because PAM finds no recurrent eligible features, the
#'   function returns `NULL`. In every successfully completed case, the
#'   complete membership is written to `subgroup_membership.csv` and saved in
#'   `subgroup_detection.rds`.
#'
#' @details Estimation is computationally intensive and writes intermediate and
#'   final artifacts to `directory`. Every successfully completed run writes
#'   `GStruc.RDS`, `subgroup_detection.rds`, `subgroup_membership.csv`, and
#'   `Subgroups Plot.png`, as well as a subgroup structure file for each
#'   membership group and a shared individual-model directory. When usable PAM
#'   distances or legacy similarities are available, the PNG is a
#'   two-dimensional subject-distance map: closer nodes have more similar
#'   individual model evidence, nearest neighbors are connected, and colors
#'   and enclosing hulls indicate membership. Subjects lacking a diagnostic
#'   distance are noted but omitted from the map. A complete membership roster
#'   is plotted when no usable pairwise distance is available.
#'
#'   With `subgroup.model = TRUE`, each subgroup directory
#'   `Models/Subgroup <g>/` additionally contains `Subgroup <g> Params.png`,
#'   one `Subgroup <g> Delta_t = <delta>.png` transition plot for every value in
#'   `time.intervals`, and `Subgroup_<g>Model.RDS`. The RDS contains one fitted
#'   OpenMx model with shared subgroup parameters and one internal likelihood
#'   block per person; those blocks are not separately estimated subject
#'   models. The transition plots display `exp(A * delta)`, the discrete-time
#'   form of the fitted continuous-time drift at the requested elapsed time.
#'
#'   PAM subgrouping requires at least three usable subject model and
#'   modification-index pairs plus at least one eligible off-diagonal,
#'   non-group path that passes its fixed recurrence screen. It chooses among
#'   candidate counts from two through `max.subgroups`, bounded above by the
#'   number of usable subjects minus one, using average silhouette width.
#'   Legacy subgrouping requires at least two usable subject model and
#'   modification-index pairs. If these requirements are not met, all subjects
#'   are assigned to membership group 1. When only some subjects have usable
#'   artifacts, either detector assigns the remaining subjects to group 1 in
#'   the complete membership output. The method-native memberships can cover
#'   only the usable subjects: PAM's `$clustering` vector and the legacy
#'   communities object omit unusable subjects, although returned PAM objects
#'   also receive the complete `$membership` vector described above.
#'
#'   `ME.free` and `PE.free` affect diagonal variances only;
#'   measurement-error and process-noise covariances remain fixed to zero.
#'   Estimating both variance sets generally requires more within-person
#'   information than estimating process noise while fixing measurement error.
#'   Estimated values are stored in the `R` and `Q` matrices, respectively, of
#'   saved OpenMx models. Internal OpenMx fits use standard console reporting
#'   rather than the interactive progress callback, which avoids
#'   `imxReportProgress` lookup failures reported in some RStudio/OpenMx
#'   installations.
#'
#' @references Park, J. J., Fisher, Z. F., Hunter, M. D., Shenk, C., Russell,
#'   M., Molenaar, P. C. M., & Chow, S.-M. (2025). Unsupervised model
#'   construction in continuous-time. *Structural Equation Modeling: A
#'   Multidisciplinary Journal, 32*(3), 377--399.
#'   \doi{10.1080/10705511.2024.2429544}
#'
#' @examples
#' quick <- ctsgimme_demo()
#' table(quick$membership)
#'
#' \dontrun{
#' fit <- ctsgimme(
#'   varnames = c("x1", "x2", "x3"),
#'   dataframe = observations,
#'   id = "id",
#'   time = "time",
#'   directory = "ctsgimme-output",
#'   sub.sig.thrsh = 0.55,
#'   max.subgroups = 4,
#'   ME.var = 0.05,
#'   ME.free = FALSE,
#'   PE.var = 1,
#'   PE.free = TRUE
#' )
#' }
#'
#' @seealso [ctsgimme_demo()] for a fast, estimation-free demonstration of
#'   the subgroup selection step.
#'
#' @export
ctsgimme <- function(varnames = NULL, dataframe = NULL,
                     id = NULL, time = NULL,
                     cores = 1, directory = NULL,
                     sig.thrsh = 0.55, sub.sig.thrsh = 1.00,
                     Galpha = 0.05, ben.hoch = TRUE, S.Galpha = 0.05, Ialpha = 0.01,
                      ME.var = 1e-8, PE.var = NULL,
                      subgroup.model = FALSE,
                      time.intervals = c(1),
                      keep.intermediate = FALSE,
                     conduct = TRUE,
                     subgroup.method = c("pam", "legacy"),
                     max.subgroups = NULL,
                      scale.data = FALSE,
                      ME.free = FALSE,
                      PE.free = FALSE) {
  subgroup.method <- match.arg(subgroup.method)
  if (length(scale.data) != 1L || !is.logical(scale.data) || is.na(scale.data)) {
    stop("scale.data must be supplied as TRUE or FALSE.")
  }
  if (length(subgroup.model) != 1L || !is.logical(subgroup.model) ||
      is.na(subgroup.model)) {
    stop("subgroup.model must be supplied as TRUE or FALSE.")
  }
  if (isTRUE(subgroup.model)) {
    .ctsgimme_validate_time_intervals(time.intervals)
  }

  .ctsgimme_validate_inputs(varnames, dataframe, id, time, directory)
  .ctsgimme_validate_subgroup_options(
    sub.sig.thrsh,
    subgroup.method,
    max.subgroups
  )
  .ctsgimme_load_dependencies(
    conduct,
    sub.sig.thrsh,
    subgroup.method,
    subgroup.model
  )
  previous_threads <- tryCatch(
    OpenMx::mxOption(NULL, "Number of Threads"),
    error = function(e) NULL
  )
  if (!is.null(previous_threads)) {
    OpenMx::mxOption(NULL, "Number of Threads", 1L)
    on.exit(
      OpenMx::mxOption(NULL, "Number of Threads", previous_threads),
      add = TRUE
    )
  }
  context <- .ctsgimme_prepare_context(
    varnames,
    dataframe,
    id,
    time,
    cores,
    directory,
    ME.var,
    PE.var,
    ben.hoch,
    ME.free,
    PE.free
  )
  if (isTRUE(scale.data)) {
    context <- .ctsgimme_scale_data(context)
  }

  worker_cluster <- .ctsgimme_make_worker_cluster(context$cores)
  if (!is.null(worker_cluster)) {
    on.exit(.ctsgimme_stop_worker_cluster(worker_cluster), add = TRUE)
  }

  G.DRIFT <- .ctsgimme_run_group_stage(
    context,
    Galpha,
    sig.thrsh,
    worker_cluster = worker_cluster
  )
  subgroup_result <- .ctsgimme_detect_subgroups(
    context,
    G.DRIFT,
    sub.sig.thrsh,
    conduct,
    subgroup.method,
    max.subgroups
  )
  memb <- subgroup_result$membership
  membership_artifacts <- .ctsgimme_write_subgroup_membership(
    memb,
    context$directory,
    distance = subgroup_result$diagnostics$distance,
    similarity = subgroup_result$diagnostics$similarity,
    diagnostic_ids = subgroup_result$diagnostics$kept.ids,
    method = subgroup_result$method
  )
  subgroup_result$diagnostics$membership.artifacts <- membership_artifacts
  .ctsgimme_save_rds(
    subgroup_result,
    file.path(context$directory, "subgroup_detection.rds")
  )
  clustering_result <- subgroup_result$clustering
  if (identical(subgroup.method, "pam") && !is.null(clustering_result)) {
    # Keep the standard PAM $clustering vector and also expose the full,
    # ID-aligned membership for direct lookup by subject identifier.
    clustering_result$membership <- memb
    attr(clustering_result, "ctsgimme.version") <- "0.0.11"
  }
  if (!is.null(clustering_result)) {
    attr(clustering_result, "ctsgimme.membership.artifacts") <-
      membership_artifacts
  }
  print(memb)

  .ctsgimme_run_subgroup_stages(
    context,
    memb,
    G.DRIFT,
    S.Galpha,
    sub.sig.thrsh,
    subgroup.model,
    time.intervals,
    Ialpha,
    worker_cluster = worker_cluster
  )
  .ctsgimme_cleanup(context, keep.intermediate)

  message(paste0(
    "Subgrouping with Continuous-Time GIMME Complete. Find networks in ",
    directory,
    "."
  ))

  if (sub.sig.thrsh == 1) {
    invisible(list(
      message = "Continuous-Time S-GIMME Complete.",
      G.DRIFT = G.DRIFT,
      membership = memb,
      directory = directory
    ))
  } else {
    clustering_result
  }
}
