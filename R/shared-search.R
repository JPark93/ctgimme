.ctgimme_fit_subjects <- function(
    context, data_for_fit, subject_ids, DRIFT, model_dir_out,
    save_models = TRUE,
    mi_dir_out = file.path(context$directory, "MIs"),
    worker_cluster = NULL) {
  # Force promises before serializing the worker closure. Without this, a
  # PSOCK worker can receive an unevaluated expression such as `sim` or
  # `context$dataframe` and try to resolve it in the worker's global session.
  force(context)
  force(data_for_fit)
  force(subject_ids)
  force(DRIFT)
  force(model_dir_out)
  force(save_models)
  force(mi_dir_out)
  force(worker_cluster)

  subject_ids <- as.character(subject_ids)
  dir.create(mi_dir_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(model_dir_out, recursive = TRUE, showWarnings = FALSE)

  fit_one <- function(i) {
    subset_dat <- data_for_fit[
      as.character(data_for_fit[[context$id]]) == as.character(i),
      ,
      drop = FALSE
    ]
    if (!nrow(subset_dat)) return(NULL)
    analysis_result <- tryCatch({
      .ctgimme_mx_try_hard(
        build_ou_model(
          subset_dat, DRIFT, context$nvar, context$varnames,
          context$ME.var, context$PE.var,
          context$ME.free, context$PE.free,
          model_name = "OUMod", lb = -10, ub = 10,
          time_col = context$time_col
        ),
        context$verbose
      )
    }, error = function(e) {
      warning("Error for subject ", i, ": ", e$message, call. = FALSE)
      NULL
    })
    if (!is.null(analysis_result)) {
      MIs <- tryCatch(OpenMx::mxMI(analysis_result, matrices = "A"), error = function(e) {
        warning("MI error for subject ", i, ": ", e$message, call. = FALSE)
        NULL
      })
      if (!is.null(MIs)) {
        .ctgimme_save_rds(MIs, file.path(mi_dir_out, paste0("MI_", i, ".RDS")))
      }
      if (isTRUE(save_models)) {
        .ctgimme_save_rds(
          analysis_result,
          file.path(model_dir_out, paste0("Model_", i, ".RDS"))
        )
      }
    }
    invisible(analysis_result)
  }

  if (context$cores <= 1 || length(subject_ids) <= 1) {
    invisible(lapply(subject_ids, fit_one))
  } else {
    owns_cluster <- is.null(worker_cluster)
    cl <- if (owns_cluster) {
      .ctgimme_make_worker_cluster(min(context$cores, length(subject_ids)))
    } else {
      worker_cluster
    }
    if (owns_cluster) {
      on.exit(.ctgimme_stop_worker_cluster(cl), add = TRUE)
    }
    clusterExport(
      cl,
      c(
        "fit_one", "build_ou_model", ".ctgimme_save_rds",
        ".ctgimme_mx_try_hard", ".ctgimme_resolve_initial_covariance",
        ".ctgimme_initial_covariance"
      ),
      envir = environment()
    )
    invisible(parLapply(cl, subject_ids, fit_one))
  }
}

.ctgimme_read_mi_stack <- function(context, valid_ids = NULL) {
  rdss <- list.files(
    file.path(context$directory, "MIs"),
    pattern = "\\.RDS$",
    full.names = TRUE
  )
  if (!is.null(valid_ids)) {
    artifact_ids <- .ctgimme_artifact_subject_id(rdss, "MI_")
    rdss <- rdss[
      !is.na(artifact_ids) & artifact_ids %in% as.character(valid_ids)
    ]
  }
  rdss <- mixedsort(rdss)
  agg <- aggregate_mis(rdss, context$param_names)
  list(files = agg$files, EPCs = agg$EPCs, rdss = rdss)
}

.ctgimme_delete_subject_files <- function(context, subject_ids, model_dir_out) {
  subject_ids <- as.character(subject_ids)
  unlink(
    file.path(context$directory, "MIs", paste0("MI_", subject_ids, ".RDS")),
    force = TRUE,
    expand = FALSE
  )
  unlink(
    file.path(model_dir_out, paste0("Model_", subject_ids, ".RDS")),
    force = TRUE,
    expand = FALSE
  )
}

.ctgimme_prune_shared <- function(
    context, DRIFT, model_dir_out, threshold, subject_ids = context$ids,
    protected_mask = NULL) {
  rdss1 <- file.path(
    model_dir_out,
    paste0("Model_", as.character(subject_ids), ".RDS")
  )
  rdss1 <- rdss1[file.exists(rdss1)]
  if (!length(rdss1)) return(list(DRIFT = DRIFT, pruned = FALSE))
  models <- list()
  for (prn in seq_along(rdss1)) {
    temp1 <- tryCatch(readRDS(rdss1[prn]), error = function(e) {
      warning("Failed to read ", prn, ": ", e$message, call. = FALSE)
      NULL
    })
    if (is.null(temp1)) next
    drifts <- tryCatch(
      subset(summary(temp1)$parameters, matrix == "A"),
      error = function(e) NULL
    )
    if (is.null(drifts) || !nrow(drifts)) next
    temp.mat1 <- matrix(NA, context$nvar, context$nvar)
    if ("name" %in% names(drifts)) {
      cells <- matrix(
        as.numeric(unlist(regmatches(drifts$name, gregexpr("\\d+", drifts$name)))),
        ncol = 2,
        byrow = TRUE
      )
      row_indices <- cells[, 1]
      col_indices <- cells[, 2]
    } else {
      row_indices <- as.numeric(drifts$row)
      col_indices <- as.numeric(drifts$col)
    }
    for (ii in seq_along(row_indices)) {
      est <- drifts[ii, "Estimate"]
      se <- drifts[ii, "Std.Error"]
      if (is.finite(est) && is.finite(se) && se > 0) {
        temp.mat1[row_indices[ii], col_indices[ii]] <- ifelse(
          abs(est) / se > qnorm(0.95), TRUE, FALSE
        )
      }
    }
    models[[length(models) + 1]] <- cbind(temp.mat1)
  }
  if (!length(models)) return(list(DRIFT = DRIFT, pruned = FALSE))
  arr <- simplify2array(models)
  true_count <- apply(
    arr,
    c(1, 2),
    function(x) sum(x == TRUE, na.rm = TRUE)
  ) / length(models)
  diag(true_count) <- 1.00
  if (!is.null(protected_mask)) true_count[protected_mask] <- 1.00
  true_count <- ifelse(true_count == 0, NA, true_count)
  true_count[DRIFT == "0"] <- NA
  if (!any(true_count <= threshold * 0.90, na.rm = TRUE)) {
    return(list(DRIFT = DRIFT, pruned = FALSE))
  }
  idx <- which(true_count == min(true_count, na.rm = TRUE))[1]
  if (!length(idx) || !is.finite(true_count[idx])) {
    return(list(DRIFT = DRIFT, pruned = FALSE))
  }
  cells <- arrayInd(idx, dim(true_count))[1, ]
  DRIFT[cells[1], cells[2]] <- "0"
  .ctgimme_inform(
    context$verbose,
    paste0("Pruned drift parameter A[", cells[1], ",", cells[2], "]")
  )
  list(DRIFT = DRIFT, pruned = TRUE)
}

.ctgimme_shared_search <- function(
    context, subject_ids, DRIFT, alpha, threshold, model_dir_out,
    subgroup = NULL, protected_mask = NULL, worker_cluster = NULL) {
  iterate <- 0
  count <- 1
  bh_thresholds <- make_threshold_sequence(
    context$nvar^2,
    alpha,
    ben_hoch = context$ben.hoch
  )

  while (iterate < 1) {
    stack <- .ctgimme_read_mi_stack(context, valid_ids = subject_ids)
    files <- stack$files
    if (is.null(files) || !ncol(files)) break

    files <- abs(files)
    usable_subjects <- colSums(is.finite(files)) > 0
    files <- files[, usable_subjects, drop = FALSE]
    n_converge <- ncol(files)
    n_requested <- length(subject_ids)
    if (is.null(files) || !n_converge) break

    if (n_converge <= (n_requested / 2)) {
      .ctgimme_inform(
        context$verbose,
        "No shared path added: usable MI files for <= half of subjects."
      )
      iterate <- 1
      next
    }

    has_any_mi <- rowSums(is.finite(files)) > 0
    if (!any(has_any_mi)) break

    # Aggregate each candidate across subjects with the median. This limits
    # the influence of unusually large, unstable modification indices.
    median_mi <- rep(-Inf, nrow(files))
    names(median_mi) <- rownames(files)
    median_mi[has_any_mi] <- matrixStats::rowMedians(
      files[has_any_mi, , drop = FALSE],
      na.rm = TRUE
    )

    best_idx <- which.max(median_mi)
    if (!length(best_idx) || !is.finite(median_mi[best_idx])) break

    selected_param <- names(median_mi)[best_idx]
    cells <- .ctgimme_get_cells(selected_param)
    if (anyNA(cells)) {
      warning(
        "Stopping shared search: malformed median-best parameter name: ",
        selected_param,
        call. = FALSE
      )
      break
    }

    n_alpha <- if (nrow(bh_thresholds) > 0) {
      bh_thresholds[min(count, nrow(bh_thresholds)), 1]
    } else {
      alpha
    }

    selected_mi <- files[selected_param, , drop = TRUE]
    selected_p <- pchisq(selected_mi, df = 1, lower.tail = FALSE)
    sig_vec <- is.finite(selected_p) & selected_p < n_alpha
    sig_count <- sum(sig_vec, na.rm = TRUE)
    sig_prop <- sig_count / n_converge
    median_p <- pchisq(median_mi[best_idx], df = 1, lower.tail = FALSE)

    if (is.finite(sig_prop) && sig_prop >= threshold) {
      DRIFT[cells[1], cells[2]] <- paste0("A_", cells[1], ",", cells[2])
      .ctgimme_inform(context$verbose, paste0(
        "Adding drift parameter A[", cells[1], ",", cells[2], "]",
        " | median MI = ", round(median_mi[best_idx], 4),
        " | significant in ", sig_count, " of ", n_converge,
        " usable subjects (", round(sig_prop, 3), ")",
        " | adapted BH alpha = ", signif(n_alpha, 4)
      ))
      .ctgimme_inform(context$verbose, paste0("Completed step ", count, "."))
      .ctgimme_delete_subject_files(context, subject_ids, model_dir_out)
      .ctgimme_fit_subjects(
        context,
        context$dataframe,
        subject_ids,
        DRIFT,
        model_dir_out,
        save_models = TRUE,
        worker_cluster = worker_cluster
      )
      count <- count + 1
    } else {
      .ctgimme_inform(context$verbose, paste0(
        "Stopping shared search: median-best candidate A[",
        cells[1], ",", cells[2], "]",
        " did not meet threshold | median MI = ", round(median_mi[best_idx], 4),
        " | median-MI p = ", signif(median_p, 4),
        " | significant in ", sig_count, " of ", n_converge,
        " usable subjects (", round(sig_prop, 3), ")",
        " | required proportion >= ", threshold,
        " | adapted BH alpha = ", signif(n_alpha, 4)
      ))
      .ctgimme_inform(context$verbose, "Beginning shared-path pruning.")
      prune <- 0
      while (prune < 1) {
        if (sum(DRIFT != "0") == context$nvar) {
          prune <- 1
          break
        }
        pruned <- .ctgimme_prune_shared(
          context,
          DRIFT,
            model_dir_out,
            threshold,
            subject_ids = subject_ids,
            protected_mask = protected_mask
        )
        DRIFT <- pruned$DRIFT
        if (!isTRUE(pruned$pruned)) {
          prune <- 1
        } else {
          .ctgimme_delete_subject_files(context, subject_ids, model_dir_out)
          .ctgimme_fit_subjects(
            context,
            context$dataframe,
            subject_ids,
            DRIFT,
            model_dir_out,
            save_models = TRUE,
            worker_cluster = worker_cluster
          )
        }
      }
      iterate <- 1
    }
  }
  DRIFT
}
