# Fit and save one subject model
fit_and_save_worker <- function(i, data, id_col, varnames, nvar, ME.var, PE.var,
                                drift_matrix, directory, lb = -Inf, ub = Inf,
                                save_model = FALSE, model_dir = NULL, scale_data = FALSE,
                                mi_extra_slash = TRUE, time_col = "Time",
                                ME.free = FALSE, PE.free = FALSE) {
  subset_dat = data[data[[id_col]] == i, , drop = FALSE]
  if (scale_data) {
    subset_dat[, varnames] = ctsgimme_safe_scale(subset_dat[, varnames, drop = FALSE])
  }
  osc = build_ou_model(subset_dat, drift_matrix, nvar, varnames, ME.var, PE.var,
                       ME.free, PE.free,
                       model_name = "OUMod", lb = lb, ub = ub, time_col = time_col)
  analysis_result = tryCatch({
    fit = OpenMx::mxTryHard(osc, silent = FALSE)
  }, error = function(e) {
    message("Error for subject ", i, ": ", e$message)
  })
  if (!is.null(analysis_result)) {
    MIs = tryCatch(OpenMx::mxMI(analysis_result, matrices = "A"), error = function(e) {
      message("Error computing MIs for subject ", i, ": ", e$message)
      NULL
    })
    if (is.null(MIs)) return(NULL)
    if (mi_extra_slash) {
      saveRDS(object = MIs, file = paste0(directory, "/MIs/", "/MI_", i, ".RDS"))
    } else {
      saveRDS(object = MIs, file = paste0(directory, "/MIs/MI_", i, ".RDS"))
    }
    if (save_model && !is.null(model_dir)) {
      saveRDS(object = analysis_result, file = paste0(model_dir, "/Model_", i, ".RDS"))
    }
  }
}

# Run subject fits with PSOCK workers
run_parallel_fit <- function(ids, data, id_col, varnames, nvar, ME.var, PE.var,
                             drift_matrix, directory, cores, time_col = "Time", lb = -Inf, ub = Inf,
                             save_model = FALSE, model_dir = NULL, scale_data = FALSE,
                             mi_extra_slash = TRUE,
                             packages = list('OpenMx', 'qgraph'),
                             ME.free = FALSE, PE.free = FALSE) {
  force(ids)
  force(data)
  force(id_col)
  force(varnames)
  force(nvar)
  force(ME.var)
  force(PE.var)
  force(ME.free)
  force(PE.free)
  force(drift_matrix)
  force(directory)
  force(cores)
  force(time_col)
  force(lb)
  force(ub)
  force(save_model)
  force(model_dir)
  force(scale_data)
  force(mi_extra_slash)
  force(packages)

  cl = makeCluster(cores, type = "PSOCK")
  on.exit(.ctsgimme_stop_worker_cluster(cl), add = TRUE)
  cl = .ctsgimme_register_worker_cluster(cl)
  clusterExport(cl, c("build_ou_model", "fit_and_save_worker", "ctsgimme_safe_scale", "packages"), envir = environment())
  clusterEvalQ(cl, {
    invisible(lapply(packages, require, character.only = TRUE))
    OpenMx::mxOption(NULL, "Number of Threads", 1L)
  })
  result = parLapply(cl, ids, function(i) {
    fit_and_save_worker(
      i, data, id_col, varnames, nvar, ME.var, PE.var,
      drift_matrix, directory, lb, ub, save_model, model_dir, scale_data,
      mi_extra_slash, time_col, ME.free, PE.free
    )
  })
  invisible(result)
}

# Aggregate modification indices
aggregate_mis <- function(rdss, param_names) {
  if (length(rdss) == 0) {
    return(list(files = NULL, EPCs = NULL))
  }
  n = length(param_names)
  nfiles = length(rdss)
  files = matrix(NA_real_, nrow = n, ncol = nfiles)
  rownames(files) = param_names
  EPCs = matrix(NA_real_, nrow = n, ncol = nfiles)
  rownames(EPCs) = param_names
  for (idx in seq_along(rdss)) {
    file = rdss[idx]
    mi_full = abs(safe_read_vector(file, "MI.Full", param_names))
    epc = safe_read_vector(file, "EPC", param_names)
    files[, idx] = mi_full
    EPCs[, idx] = epc
  }
  return(list(files = files, EPCs = EPCs))
}

# Compute pruning proportions
compute_pruning_proportions <- function(model_dir, nvar) {
  rdss1 = list.files(model_dir, pattern = "\\.RDS$", full.names = TRUE)
  if (length(rdss1) == 0) {
    return(NULL)
  }
  models = list()
  for (prn in seq_along(rdss1)) {
    temp1 = tryCatch({readRDS(rdss1[prn])}, error = function(e) {
      message("Failed to read ", prn, ": ", e$message)
      NULL})
    if (is.null(temp1)) next
    params = tryCatch(summary(temp1)$parameters, error = function(e) {
      message("Failed to summarize ", prn, ": ", e$message)
      NULL
    })
    if (is.null(params) || nrow(params) == 0) next
    drifts = subset(params, matrix == 'A')
    if (nrow(drifts) == 0) next
    cells = matrix(
      as.numeric(unlist(regmatches(drifts$name, gregexpr("\\d+", drifts$name)))),
      ncol = 2, byrow = TRUE
    )
    row_indices = cells[, 1]
    col_indices = cells[, 2]
    temp.mat1 = matrix(NA, nvar, nvar)
    for (i in 1:nrow(cells)) {
      temp.mat1[row_indices[i], col_indices[i]] = ifelse(abs(drifts[i, "Estimate"]) /
                                                           (drifts[i, "Std.Error"]) >
                                                           qnorm(0.975), TRUE, FALSE)
    }
    models[[length(models) + 1]] = cbind(temp.mat1)
  }
  if (length(models) == 0) {
    return(NULL)
  }
  arr = simplify2array(models)
  true_count = apply(arr, c(1, 2), function(x) sum(x == TRUE, na.rm = TRUE)) / length(models)
  return(true_count)
}

# Plot group paths
plot_group_paths <- function(G.DRIFT, varnames, directory) {
  output_path = file.path(directory, "Group Paths.png")
  png(filename = output_path, width = 800, height = 800)
  qgraph(t(G.DRIFT != "0") * 1, layout = "circle", labels = varnames,
         edge.width = 5, diag = TRUE, edge.labels = "GROUP")
  dev.off()
}
