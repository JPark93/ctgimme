.ctsgimme_fit_individual <- function(
    context, i, new.data, SG.DRIFT, m, nks, Ialpha) {
  DRIFT <- SG.DRIFT
  subset_dat <- new.data[
    as.character(new.data[[context$id]]) == as.character(i),
    ,
    drop = FALSE
  ]
  osc <- build_ou_model(
    subset_dat,
    DRIFT,
    context$nvar,
    context$varnames,
    context$ME.var,
    context$PE.var,
    context$ME.free,
    context$PE.free,
    model_name = "OUMod",
    lb = -10,
    ub = 10,
    time_col = context$time_col
  )
  fit <- tryCatch(OpenMx::mxTryHard(osc, silent = FALSE), error = function(e) {
    message("Error for subject ", i, ": ", e$message)
    NULL
  })
  if (is.null(fit)) return(NULL)
  fit2 <- fit
  optimization <- 0
  count <- 1
  while (optimization < 1) {
    if (count > 1) {
      fit2 <- fit
      fit <- tryCatch(
        OpenMx::mxTryHard(osc, silent = FALSE),
        error = function(e) NULL
      )
      if (is.null(fit) || any(is.na(summary(fit)$parameters$Std.Error))) {
        fit <- fit2
        optimization <- 1
        break
      }
    }
    if (m <= 0 || nrow(nks) <= 0 || count > nrow(nks)) {
      optimization <- 1
      break
    }
    MIs <- tryCatch(OpenMx::mxMI(fit, matrices = "A"), error = function(e) NULL)
    if (is.null(MIs) || is.null(MIs$MI.Full)) {
      optimization <- 1
      fit <- fit2
    } else if (length(MIs$MI.Full) == length(names(MIs$MI.Full))) {
      max_idx <- which.max(abs(MIs$MI.Full))
      max_val <- abs(MIs$MI.Full[max_idx])
      if (is.finite(max_val) &&
          max_val >= qchisq(1 - nks[count, ], df = 1)) {
        max_name <- names(MIs$MI.Full)[max_idx]
        cells <- .ctsgimme_get_cells(max_name)
        osc$A$free[cells[1], cells[2]] <- TRUE
        osc$A$labels[cells[1], cells[2]] <- NA_character_
        message(paste0(
          "Adding drift parameter A[", cells[1], ",", cells[2], "]"
        ))
        count <- count + 1
        if (sum(osc$A$free) == context$nvar^2) optimization <- 1
      } else {
        optimization <- 1
      }
    } else {
      message(paste0(
        "Malformed MI.Full: names and values mismatch. Skipping subject ",
        i
      ))
      optimization <- 1
      fit <- fit2
    }
  }

  message("Pruning Stage.")
  prune <- 0
  prune_iter <- 0
  max_prune_iter <- context$nvar^2
  while (prune < 1) {
    prune_iter <- prune_iter + 1
    if (prune_iter > max_prune_iter) {
      message("Pruning stopped after maximum number of iterations.")
      break
    }
    stat.sig <- tryCatch(summary(fit)$parameters, error = function(e) NULL)
    if (is.null(stat.sig) || !nrow(stat.sig)) break
    stat.sig <- subset(stat.sig, matrix == "A")
    if (!nrow(stat.sig)) break
    stat.sig$z <- abs(stat.sig$Estimate) / stat.sig$Std.Error
    idx <- cbind(as.integer(stat.sig$row), as.integer(stat.sig$col))
    protected <- SG.DRIFT != "0"
    prunable <- stat.sig[!protected[idx], , drop = FALSE]
    prunable <- prunable[is.finite(prunable$z), , drop = FALSE]
    if (nrow(prunable) == 0 || is.null(nrow(prunable))) {
      message("No unprotected parameters left to prune.")
      break
    }
    min_z_index <- which.min(prunable$z)
    if (length(prunable$z[min_z_index]) == 0) {
      message("All Prunable Paths Removed.")
      break
    }
    if (prunable$z[min_z_index] < qnorm(1 - Ialpha)) {
      this <- prunable[min_z_index, ]
      cells <- matrix(c(this$row, this$col), 1, 2)
      osc$A$free[cells[1, 1], cells[1, 2]] <- FALSE
      osc$A$labels[cells[1, 1], cells[1, 2]] <- NA_character_
      new_fit <- tryCatch(
        OpenMx::mxTryHard(osc, silent = FALSE),
        error = function(e) NULL
      )
      if (is.null(new_fit)) {
        osc$A$free[cells[1, 1], cells[1, 2]] <- TRUE
        message("Pruning fit failed; retaining last converged fit.")
        break
      }
      fit <- new_fit
      message(paste0(
        "NOTE: Pruned drift parameter: A[",
        cells[1, 1], ",", cells[1, 2], "]!"
      ))
    } else {
      message("No Pruning Conducted.")
      prune <- 1
    }
  }

  if (!is.null(fit)) {
    .ctsgimme_save_rds(
      fit,
      file.path(
        context$directory,
        "Models",
        "Individuals",
        paste0("FinalModel_", i, ".RDS")
      )
    )
    stat.sig <- subset(summary(fit)$parameters, matrix == "A")
    ests <- matrix(0, context$nvar, context$nvar)
    for (jj in seq_len(nrow(stat.sig))) {
      ests[stat.sig$row[jj], stat.sig$col[jj]] <- stat.sig$Estimate[jj]
    }
    ests <- t(ests)
    .ctsgimme_safe_png(
      file.path(
        context$directory,
        "Models",
        "Individuals",
        paste0("FinalModel_", i, ".PNG")
      ),
      {
        qgraph(
          ests,
          layout = "circle",
          labels = context$varnames,
          edge.width = 1,
          diag = TRUE,
          edge.labels = round(ests, 2),
          theme = "colorblind",
          fade = FALSE
        )
      }
    )
  }
  invisible(fit)
}

.ctsgimme_run_individual_fits <- function(
    context, valid_ids, new.data, SG.DRIFT, m, nks, Ialpha,
    worker_cluster = NULL) {
  # Resolve caller promises on the master before exporting the closure to
  # PSOCK workers; worker sessions cannot resolve caller-local symbols.
  force(context)
  force(valid_ids)
  force(new.data)
  force(SG.DRIFT)
  force(m)
  force(nks)
  force(Ialpha)
  force(worker_cluster)

  fit_individual <- function(i) {
    .ctsgimme_fit_individual(
      context,
      i,
      new.data,
      SG.DRIFT,
      m,
      nks,
      Ialpha
    )
  }

  if (context$cores <= 1 || length(valid_ids) <= 1) {
    invisible(lapply(valid_ids, fit_individual))
  } else {
    owns_cluster <- is.null(worker_cluster)
    cl <- if (owns_cluster) {
      .ctsgimme_make_worker_cluster(min(context$cores, length(valid_ids)))
    } else {
      worker_cluster
    }
    if (owns_cluster) {
      on.exit(.ctsgimme_stop_worker_cluster(cl), add = TRUE)
    }
    clusterExport(
      cl,
      c(
        "fit_individual", ".ctsgimme_fit_individual",
        "build_ou_model", ".ctsgimme_save_rds", ".ctsgimme_safe_png",
        ".ctsgimme_get_cells"
      ),
      envir = environment()
    )
    invisible(parLapply(cl, valid_ids, fit_individual))
  }
}

.ctsgimme_run_subgroup_stages <- function(
    context, memb, G.DRIFT, S.Galpha, sub.sig.thrsh, subgroup.model,
    time.intervals, Ialpha, worker_cluster = NULL) {
  for (subgroup in sort(unique(memb))) {
    DRIFT <- G.DRIFT
    saveRDS(DRIFT, file.path(context$directory, "GStruc.RDS"))
    subgroup_dir <- file.path(
      context$directory,
      "Models",
      paste0("Subgroup ", subgroup)
    )
    dir.create(subgroup_dir, showWarnings = FALSE, recursive = TRUE)
    sg.ids <- names(memb)[memb == subgroup]
    new.data <- context$dataframe[
      as.character(context$dataframe[[context$id]]) %in% as.character(sg.ids),
      ,
      drop = FALSE
    ]
    valid_ids <- unique(new.data[[context$id]])
    if (!length(valid_ids)) next

    DRIFT <- .ctsgimme_shared_search(
      context,
      valid_ids,
      DRIFT,
      S.Galpha,
      sub.sig.thrsh,
      subgroup_dir,
      subgroup = subgroup,
      protected_mask = G.DRIFT != "0",
      worker_cluster = worker_cluster
    )
    message(paste0(
      "Subgroup Search ", subgroup, " of ", max(memb), " Complete."
    ))

    if (subgroup.model == TRUE) {
      .ctsgimme_fit_subgroup_model(
        context,
        new.data,
        DRIFT,
        G.DRIFT,
        subgroup,
        subgroup_dir,
        time.intervals
      )
    }

    message("Beginning Individual Model Fitting for Subgroup Members")
    m <- (context$nvar^2) - sum(DRIFT != "0")
    nks <- make_threshold_sequence(
      context$nvar^2,
      Ialpha,
      ben_hoch = context$ben.hoch
    )
    n_higher_paths <- sum(DRIFT != "0")
    if (nrow(nks) > 0 && n_higher_paths > 0) {
      if (n_higher_paths >= nrow(nks)) {
        nks <- matrix(numeric(0), 0, 1)
      } else {
        nks <- nks[-seq_len(n_higher_paths), , drop = FALSE]
      }
    }
    SG.DRIFT <- DRIFT
    saveRDS(
      SG.DRIFT,
      file.path(context$directory, paste0("SGStruc", subgroup, ".RDS"))
    )
    .ctsgimme_safe_png(
      file.path(
        subgroup_dir,
        paste0("Subgroup ", subgroup, " Paths.png")
      ),
      {
        qgraph(
          t(abs(((G.DRIFT != "0") * 1) - ((DRIFT != "0") * 1))),
          layout = "circle",
          labels = context$varnames,
          edge.width = 5,
          diag = TRUE,
          edge.labels = paste0("SG-", subgroup)
        )
      }
    )

    .ctsgimme_run_individual_fits(
      context,
      valid_ids,
      new.data,
      SG.DRIFT,
      m,
      nks,
      Ialpha,
      worker_cluster = worker_cluster
    )
  }
  invisible(NULL)
}
