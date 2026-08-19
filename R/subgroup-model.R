.ctsgimme_detect_subgroups <- function(
    context, group_drift, sub.sig.thrsh, conduct,
    subgroup.method = c("pam", "legacy"), max.subgroups = NULL) {
  subgroup.method <- match.arg(subgroup.method)
  .ctsgimme_validate_subgroup_options(
    sub.sig.thrsh,
    subgroup.method,
    max.subgroups
  )
  if (sub.sig.thrsh == 1.00) {
    return(.ctsgimme_one_subgroup_result(
      context,
      method = subgroup.method,
      reason = "Subgroup detection was disabled because sub.sig.thrsh = 1."
    ))
  }
  if (identical(subgroup.method, "pam")) {
    return(.ctsgimme_detect_subgroups_pam(
      context,
      group_drift = group_drift,
      max.subgroups = max.subgroups
    ))
  }
  .ctsgimme_detect_subgroups_legacy(
    context,
    sub.sig.thrsh = sub.sig.thrsh,
    conduct = conduct
  )
}

.ctsgimme_edge_label_matrix <- function(edge_labels, nvar) {
  if (length(edge_labels) != nvar^2) {
    stop("edge_labels must contain exactly nvar^2 values.")
  }
  matrix(edge_labels, nrow = nvar, ncol = nvar)
}

.ctsgimme_subgroup_plot_components <- function(params.A, DRIFT, G.DRIFT, nvar) {
  effects <- matrix(0, nvar, nvar)
  edge_labels <- matrix("", nvar, nvar)
  significance <- ifelse(
    abs(params.A$Estimate) > qnorm(0.95) * params.A$Std.Error,
    "*",
    "ns"
  )
  for (index in seq_len(nrow(params.A))) {
    plot_row <- params.A$col[[index]]
    plot_col <- params.A$row[[index]]
    effects[plot_row, plot_col] <- params.A$Estimate[[index]]
    edge_labels[plot_row, plot_col] <- paste0(
      significance[[index]],
      " (",
      round(params.A$Estimate[[index]], 2),
      ")"
    )
  }

  shared <- t(G.DRIFT != "0" & DRIFT != "0")
  subgroup_only <- t(G.DRIFT == "0" & DRIFT != "0")
  absent <- t(DRIFT == "0")
  edge_colors <- matrix("black", nvar, nvar)
  edge_colors[shared] <- "gray"
  edge_colors[subgroup_only & effects > 0] <- "blue"
  edge_colors[subgroup_only & effects < 0] <- "red"
  edge_colors[subgroup_only & effects == 0] <- "black"
  edge_colors[absent] <- "black"

  list(
    effects = effects,
    edge.labels = edge_labels,
    edge.colors = c(edge_colors)
  )
}

.ctsgimme_fit_subgroup_model <- function(
    context, new.data, DRIFT, G.DRIFT, subgroup, subgroup_dir,
    time.intervals) {
  .ctsgimme_validate_time_intervals(time.intervals)
  message(paste0(
    "Fitting Parameterized Model of Subgroup ",
    subgroup,
    " with one multisubject likelihood"
  ))
  fit <- tryCatch(
    {
      model <- build_multisubject_ou_model(
        data = new.data,
        id_col = context$id,
        drift_matrix = DRIFT,
        nvar = context$nvar,
        varnames = context$varnames,
        ME.var = context$ME.var,
        PE.var = context$PE.var,
        ME.free = context$ME.free,
        PE.free = context$PE.free,
        model_name = "OUMod",
        lb = -10,
        ub = 10,
        time_col = context$time_col
      )
      OpenMx::mxTryHard(
        model,
        silent = FALSE
      )
    },
    error = function(e) {
      message(
        "Subgroup parameterized model failed for subgroup ",
        subgroup,
        ": ",
        e$message
      )
      NULL
    }
  )
  if (is.null(fit)) return(invisible(NULL))

  .ctsgimme_save_rds(
    fit,
    file.path(subgroup_dir, paste0("Subgroup_", subgroup, "Model.RDS"))
  )

  sum.fit <- summary(fit)
  params.A <- subset(sum.fit$parameters, matrix == "A")
  plot_components <- .ctsgimme_subgroup_plot_components(
    params.A,
    DRIFT,
    G.DRIFT,
    context$nvar
  )
  effects <- plot_components$effects
  .ctsgimme_safe_png(
    file.path(
      subgroup_dir,
      paste0("Subgroup ", subgroup, " Params.png")
    ),
    {
      qgraph(
        effects,
        layout = "circle",
        labels = context$varnames,
        edge.width = 1,
        diag = TRUE,
        edge.labels = plot_components$edge.labels,
        theme = "colorblind",
        edge.color = plot_components$edge.colors,
        fade = TRUE,
        mar = rep(8, 4)
      )
    }
  )
  if (!requireNamespace("expm", quietly = TRUE)) {
    stop(
      "Package 'expm' is required for subgroup delta-t plots when ",
      "subgroup.model = TRUE."
    )
  }
  for (ints in time.intervals) {
    delt <- round(expm::expm(effects * ints), 3)
    .ctsgimme_safe_png(
      file.path(
        subgroup_dir,
        paste0("Subgroup ", subgroup, " Delta_t = ", ints, ".png")
      ),
      {
        qgraph(
          delt,
          layout = "circle",
          labels = context$varnames,
          fade = TRUE,
          edge.width = 1,
          diag = TRUE,
          edge.labels = delt,
          maximum = 1.00,
          theme = "colorblind",
          title = paste0("Subgroup ", subgroup, "; Delta_t = ", ints),
          mar = c(8, 8, 10, 8)
        )
      }
    )
  }
  invisible(fit)
}
