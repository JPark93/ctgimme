.ctsgimme_scale_data <- function(context) {
  message("Scaling Variables for Analysis")
  for (ID in context$ids) {
    idx <- which(
      as.character(context$dataframe[[context$id]]) == as.character(ID)
    )
    context$dataframe[idx, context$varnames] <- ctsgimme_safe_scale(
      context$dataframe[idx, context$varnames, drop = FALSE]
    )
  }
  context
}

.ctsgimme_initial_drift <- function(nvar) {
  DRIFT <- matrix("0", nrow = nvar, ncol = nvar)
  diag(DRIFT) <- paste0("A_", seq_len(nvar), ",", seq_len(nvar))
  DRIFT
}

.ctsgimme_run_group_stage <- function(
    context, Galpha, sig.thrsh, worker_cluster = NULL) {
  model_dir <- file.path(context$directory, "Models")
  .ctsgimme_delete_subject_files(
    context,
    context$ids,
    model_dir
  )

  DRIFT <- .ctsgimme_initial_drift(context$nvar)
  .ctsgimme_fit_subjects(
    context,
    context$dataframe,
    context$ids,
    DRIFT,
    model_dir,
    save_models = TRUE,
    worker_cluster = worker_cluster
  )
  DRIFT <- .ctsgimme_shared_search(
    context,
    context$ids,
    DRIFT,
    Galpha,
    sig.thrsh,
    model_dir,
    worker_cluster = worker_cluster
  )
  G.DRIFT <- DRIFT
  saveRDS(G.DRIFT, file.path(context$directory, "GStruc.RDS"))
  .ctsgimme_safe_png(file.path(context$directory, "Group Paths.png"), {
    qgraph(
      t(G.DRIFT != "0") * 1,
      layout = "circle",
      labels = context$varnames,
      edge.width = 5,
      diag = TRUE,
      edge.labels = "GROUP"
    )
  })
  message("Group Search Complete.")
  G.DRIFT
}

.ctsgimme_one_subgroup_result <- function(
    context, method, reason = NULL, warn = FALSE) {
  memb <- rep(1L, length(context$ids))
  names(memb) <- as.character(context$ids)
  if (!is.null(reason)) {
    if (isTRUE(warn)) warning(reason, call. = FALSE) else message(reason)
  }
  list(
    membership = memb,
    clustering = NULL,
    walktrap = NULL,
    method = method,
    diagnostics = list(reason = reason)
  )
}
