# ============================================================================
# ctsgimme-master.R
# Refactored from ctsgimme-ver-0.0.3.R
# Corrected master version; v0.0.5 restores original adapted-BH and mean-MI shared-path selection.
# ============================================================================

# ============================================================
# Helper Function A: JPmx (corrected modification-index helper)
# ============================================================
JPmx = function (model, matrices = NA, full = TRUE){
  OpenMx:::warnModelCreatedByOldVersion(model)
  if (OpenMx:::single.na(matrices)) {
    matrices = names(model$matrices)
    if (is(model$expectation, "MxExpectationRAM")) {
      matrices = setdiff(matrices, model$expectation$F)
    }
  }
  if (imxHasWLS(model)) {
    stop("modification indices not implemented for WLS fitfunction")
  }
  param = omxGetParameters(model)
  param.names = names(param)
  gmodel = omxSetParameters(model, free = FALSE, labels = param.names)
  mi.r = NULL
  mi.f = NULL
  epc.f = NULL
  a.names = NULL
  new.models = list()
  for (amat in matrices) {
    matObj = model[[amat]]
    freemat = matObj$free
    sym.sel = upper.tri(freemat, diag = TRUE)
    notSymDiag = !(is(gmodel[[amat]])[1] %in% c("DiagMatrix", 
                                                "SymmMatrix"))
    for (i in 1:length(freemat)) {
      if (freemat[i] == FALSE && (notSymDiag || sym.sel[i] == 
                                  TRUE)) {
        tmpLab = gmodel[[amat]]$labels[i]
        plusOneParamModel = model
        if (length(tmpLab) > 0 && !is.na(tmpLab)) {
          gmodel = omxSetParameters(gmodel, labels = tmpLab, 
                                    free = TRUE)
          plusOneParamModel = omxSetParameters(plusOneParamModel, 
                                               labels = tmpLab, free = TRUE)
        }
        else {
          gmodel[[amat]]$free[i] = TRUE
          plusOneParamModel[[amat]]$free[i] = TRUE
        }
        if (is(gmodel[[amat]])[1] %in% c("ZeroMatrix")) {
          cop = gmodel[[amat]]
          newSingleParamMat = mxMatrix("Full", nrow = nrow(cop), 
                                       ncol = ncol(cop), values = cop$values, free = cop$free, 
                                       labels = cop$labels, name = cop$name, lbound = cop$lbound, 
                                       ubound = cop$ubound, dimnames = dimnames(cop))
          bop = plusOneParamModel[[amat]]
          newPlusOneParamMat = mxMatrix("Full", nrow = nrow(bop), 
                                        ncol = ncol(bop), values = bop$values, free = bop$free, 
                                        labels = bop$labels, name = bop$name, lbound = bop$lbound, 
                                        ubound = bop$ubound, dimnames = dimnames(bop))
        }
        else if (is(gmodel[[amat]])[1] %in% c("DiagMatrix", 
                                              "SymmMatrix")) {
          cop = gmodel[[amat]]
          newSingleParamMat = mxMatrix("Symm", nrow = nrow(cop), 
                                       ncol = ncol(cop), values = cop$values, free = (cop$free | 
                                                                                        t(cop$free)), labels = cop$labels, name = cop$name, 
                                       lbound = cop$lbound, ubound = cop$ubound, 
                                       dimnames = dimnames(cop))
          bop = plusOneParamModel[[amat]]
          newPlusOneParamMat = mxMatrix("Symm", nrow = nrow(bop), 
                                        ncol = ncol(bop), values = bop$values, free = (bop$free | 
                                                                                         t(bop$free)), labels = bop$labels, name = bop$name, 
                                        lbound = bop$lbound, ubound = bop$ubound, 
                                        dimnames = dimnames(bop))
        }
        else {
          newSingleParamMat = gmodel[[amat]]
          newPlusOneParamMat = plusOneParamModel[[amat]]
        }
        gmodel[[amat]] = newSingleParamMat
        plusOneParamModel[[amat]] = newPlusOneParamMat
        custom.compute = mxComputeSequence(list(mxComputeNumericDeriv(checkGradient = FALSE), 
                                                mxComputeReportDeriv()))
        gmodel = mxModel(gmodel, custom.compute)
        grun = try(mxRun(gmodel, silent = TRUE, suppressWarnings = FALSE, 
                         unsafe = TRUE))
        if (is(grun, "try-error")) {
          gmodel = omxSetParameters(gmodel, labels = names(omxGetParameters(gmodel)), 
                                    free = FALSE)
          next
        }
        grad = grun$output$gradient
        hess = grun$output$hessian
        modind = 0.5 * grad^2/hess
        # if (full == TRUE) {
        #   custom.compute.smart = mxComputeSequence(list(mxComputeNumericDeriv(knownHessian = model$output$hessian, 
        #                                                                       checkGradient = FALSE), mxComputeReportDeriv()))
        #   plusOneParamRun = mxRun(mxModel(plusOneParamModel, 
        #                                   custom.compute.smart), silent = TRUE, suppressWarnings = FALSE, 
        #                           unsafe = TRUE)
        #   grad.full = plusOneParamRun$output$gradient
        #   grad.full[is.na(grad.full)] = 0
        #   
        #   tol = max(1e-8, 1e-6 * max(abs(grad.full), na.rm = TRUE))
        #   nz  = which(abs(grad.full) > tol)
        #   
        #   hess.full = plusOneParamRun$output$hessian
        #   modind.full = 0.5 * t(matrix(grad.full)) %*% 
        #     solve(hess.full) %*% matrix(grad.full)
        #   if (sum(grad.full != 0) == 1) {
        #     exppar.full = -modind.full/grad.full[grad.full != 0]
        #   }
        #   else {
        #     stop("Something strange in the neighborhood.\nFound a one-parameter model with more than one parameter.\nPost this to the OpenMx forums.")
        #   }
        # }
        # else {
        #   modind.full = NULL
        #   exppar.full = NULL
        # }
        # --- replace your current full==TRUE block with this ---
        if (full == TRUE) {
          custom.compute.smart = mxComputeSequence(list(
            mxComputeNumericDeriv(knownHessian = model$output$hessian, checkGradient = FALSE),
            mxComputeReportDeriv()
          ))
          plusOneParamRun = mxRun(mxModel(plusOneParamModel, custom.compute.smart),
                                  silent = TRUE, suppressWarnings = FALSE, unsafe = TRUE)
          
          grad.full = plusOneParamRun$output$gradient
          if (is.null(names(grad.full))) {
            # fall back to parameter names in the plus-one model
            names(grad.full) = names(omxGetParameters(plusOneParamRun))
          }
          grad.full[!is.finite(grad.full)] = 0
          
          hess.full = plusOneParamRun$output$hessian
          modind.full = as.numeric(0.5 * t(as.matrix(grad.full)) %*% solve(hess.full) %*% as.matrix(grad.full))
          
          # name of the single freed parameter in the single-parameter run
          n.names = names(omxGetParameters(grun))
          
          # pick the gradient entry for that parameter (with tolerance fallback)
          tol = max(1e-8, 1e-6 * max(abs(grad.full), na.rm = TRUE))
          idx = match(n.names, names(grad.full))
          if (is.na(idx)) {
            nz = which(abs(grad.full) > tol)
            idx = if (length(nz)) nz[ which.max(abs(grad.full[nz])) ] else NA_integer_
          }
          
          if (is.na(idx) || abs(grad.full[idx]) <= tol) {
            exppar.full = NA_real_
          } else {
            exppar.full = -modind.full / grad.full[idx]
          }
        } else {
          modind.full = NULL
          exppar.full = NULL
        }
        
        n.names = names(omxGetParameters(grun))
        if (length(modind) > 0) {
          a.names = c(a.names, n.names)
          mi.r = c(mi.r, modind)
          mi.f = c(mi.f, modind.full)
          epc.f = c(epc.f, exppar.full)
          new.models = c(new.models, plusOneParamModel)
        }
        gmodel = omxSetParameters(gmodel, labels = names(omxGetParameters(gmodel)), 
                                  free = FALSE)
      }
    }
    names(mi.r) = a.names
    if (full == TRUE) {
      names(mi.f) = a.names
      names(epc.f) = a.names
    }
    names(new.models) = a.names
  }
  if (length(model$submodels) > 0) {
    for (asubmodel in names(model$submodels)) {
      ret = c(ret, mxMI(asubmodel))
    }
  }
  if (is.null(mi.f)) 
    mi.f = rep(NA, length(mi.r))
  if (is.null(epc.f)) 
    epc.f = rep(NA, length(mi.r))
  retList = list2DF(list(MI = mi.r, MI.Full = mi.f, plusOneParamModels = new.models, 
                         EPC = epc.f), nrow = length(mi.r))
  return(retList)
}

# ============================================================
# Helper Function B: safe_read_vector
# ============================================================
safe_read_vector <- function(file, element, param_names) {
  out = rep(NA_real_, length(param_names))
  names(out) = param_names
  
  vec = tryCatch({
    raw = readRDS(file)[[element]]
    vec = c(raw)
    if (!is.null(names(vec)) && length(vec) != length(names(vec))) {
      stop("Length of vector and names do not match")
    }
    vec
  }, error = function(e) {
    message("Failed to read ", element, " from ", file, ": ", e$message)
    return(NULL)
  })
  
  if (!is.null(vec)) {
    intersect_names = intersect(names(vec), param_names)
    out[intersect_names] = vec[intersect_names]
  }
  
  return(out)
}


# ============================================================
# Helper Function B2: robust subject-level scaling
# ============================================================
ctsgimme_safe_scale <- function(x) {
  mat = as.matrix(x)
  means = colMeans(mat, na.rm = TRUE)
  sds = apply(mat, 2, sd, na.rm = TRUE)
  centered = sweep(mat, 2, means, "-")
  out = centered
  colnames(out) = colnames(mat)
  good = is.finite(sds) & sds > 0
  if (any(good)) {
    out[, good] = sweep(centered[, good, drop = FALSE], 2, sds[good], "/")
  }
  if (any(!good)) {
    out[, !good] = centered[, !good, drop = FALSE]
  }
  out
}

# ============================================================
# Helper Function C: compute_param_names
# ============================================================
compute_param_names <- function(nvar, model_name = "OUMod") {
  param_names = character(0)
  for (j in 1:nvar) {
    for (i in 1:nvar) {
      param_names = c(param_names, sprintf("%s.A[%d,%d]", model_name, i, j))
    }
  }
  return(param_names)
}

# ============================================================
# Helper Function D: extract_matrix_indices
# ============================================================
extract_matrix_indices <- function(names_vec) {
  nums = as.numeric(unlist(regmatches(names_vec, gregexpr("\\d+", names_vec))))
  matrix(nums, ncol = 2, byrow = TRUE)
}

# ============================================================
# Helper Function E: make_threshold_sequence
# ============================================================
make_threshold_sequence <- function(m, alpha, ben_hoch = TRUE) {
  # v0.0.5: restore the original ctgimme adapted Benjamini-Hochberg
  # construction from ctgimme2282024.R:
  #   ks[k] = (k / m) * alpha, sorted decreasing.
  # This replaces the overly strict alpha / n_converge shared-path cutoff
  # used in v0.0.4.
  if (m <= 0) {
    return(matrix(numeric(0), 0, 1))
  }
  ks = matrix(NA, m, 1)
  for (k in 1:m) {
    ks[k, ] = (k / m) * alpha
  }
  ks = matrix(sort(ks, TRUE), m, 1)
  if (!ben_hoch) {
    ks = matrix(alpha, nrow(ks), 1)
  }
  return(ks)
}

# ============================================================
# Helper Function F: build_ou_model
# ============================================================
build_ou_model <- function(data, drift_matrix, nvar, varnames, ME.var, PE.var, 
                           model_name = "OUMod", lb = -Inf, ub = Inf, time_col = "Time") {
  amat = mxMatrix("Full", nvar, nvar,
                  free   = drift_matrix != "0",
                  name   = "A",
                  ubound = ub, lbound = lb)
  bmat = mxMatrix('Zero', nvar, nvar, name='B')
  cdim = list(varnames, paste0('F', 1:nvar))
  cmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name = 'C', dimnames = cdim)
  dmat = mxMatrix('Zero', nvar, nvar, name='D')
  qmat = mxMatrix('Diag', nvar, nvar, diag(FALSE, nvar), PE.var, name='Q', ubound = ub, lbound = 1e-5)
  rmat = mxMatrix('Diag', nvar, nvar, diag(FALSE, nvar), ME.var, name='R')
  xmat = mxMatrix('Full', nvar, 1, FALSE, rep(0, nvar), name='x0')
  pmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name='P0')
  umat = mxMatrix('Zero', nvar, 1, name='u')
  tmat = mxMatrix('Full', 1, 1, FALSE, name='time', labels=paste0('data.', time_col))
  mxModel(model_name, 
          amat, bmat, cmat, dmat, qmat, 
          rmat, xmat, pmat, umat, tmat,
          mxExpectationSSCT('A', 'B', 'C', 'D', 'Q', 
                            'R', 'x0', 'P0', 'u', 'time'),
          mxFitFunctionML(),
          mxData(data, 'raw'))
}

# ============================================================
# Helper Function G: fit_and_save_worker
# ============================================================
fit_and_save_worker <- function(i, data, id_col, varnames, nvar, ME.var, PE.var, 
                                drift_matrix, directory, lb = -Inf, ub = Inf, 
                                save_model = FALSE, model_dir = NULL, scale_data = FALSE,
                                mi_extra_slash = TRUE, time_col = "Time") {
  subset_dat = data[data[[id_col]] == i, , drop = FALSE]
  if (scale_data) {
    subset_dat[, varnames] = ctsgimme_safe_scale(subset_dat[, varnames, drop = FALSE])
  }
  osc = build_ou_model(subset_dat, drift_matrix, nvar, varnames, ME.var, PE.var, 
                       model_name = "OUMod", lb = lb, ub = ub, time_col = time_col)
  analysis_result = tryCatch({
    fit = mxTryHard(osc)
  }, error = function(e) {
    message("Error for subject ", i, ": ", e$message)
  })
  if (!is.null(analysis_result)) {
    MIs = tryCatch(JPmx(analysis_result, matrices = "A"), error = function(e) {
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

# ============================================================
# Helper Function H: run_parallel_fit
# ============================================================
run_parallel_fit <- function(ids, data, id_col, varnames, nvar, ME.var, PE.var, 
                             drift_matrix, directory, cores, time_col = "Time", lb = -Inf, ub = Inf, 
                             save_model = FALSE, model_dir = NULL, scale_data = FALSE,
                             mi_extra_slash = TRUE,
                             packages = list('dynr', 'OpenMx', 'qgraph')) {
  cl = makeCluster(cores, type = "PSOCK")
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, c("JPmx", "build_ou_model", "fit_and_save_worker", "ctsgimme_safe_scale", "packages"), envir = environment())
  clusterEvalQ(cl, {
    invisible(lapply(packages, require, character.only = T))
  })
  result = parLapply(cl, ids, function(i) {
    fit_and_save_worker(i, data, id_col, varnames, nvar, ME.var, PE.var, 
                        drift_matrix, directory, lb, ub, save_model, model_dir, scale_data, mi_extra_slash, time_col)
  })
  invisible(result)
}

# ============================================================
# Helper Function I: aggregate_mis
# ============================================================
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

# ============================================================
# Helper Function J: compute_pruning_proportions
# ============================================================
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

# ============================================================
# Helper Function K: plot_group_paths
# ============================================================
plot_group_paths <- function(G.DRIFT, varnames, directory) {
  output_path = file.path(directory, "Group Paths.png")
  png(filename = output_path, width = 800, height = 800)
  qgraph(t(G.DRIFT != "0") * 1, layout = "circle", labels = varnames, 
         edge.width = 5, diag = TRUE, edge.labels = "GROUP")
  dev.off()
}

# ============================================================
# Main Function: ctsgimme
# Unified master version v0.0.5: keeps the corrected master helpers and
# ctsgimme ver 0.0.4 infrastructure, while restoring original adapted-BH
# and mean-MI shared-path selection from ctgimme2282024.R.
# ============================================================
ctsgimme <- function(varnames = NULL, dataframe = NULL,
                     id = NULL, time = NULL,
                     cores = 1, directory = NULL,
                     sig.thrsh = 0.55, sub.sig.thrsh = 1.00,
                     Galpha = 0.05, ben.hoch = TRUE, S.Galpha = 0.05, Ialpha = 0.01,
                     ME.var = 1e-8, PE.var = NULL,
                     subgroup.model = FALSE,
                     time.intervals = c(1),
                     subgroup.time.mode = c("legacy", "observed", "schedule"),
                     measurement.schedule = seq(0, by = 1/8, length.out = 5),
                     cycle.interval = 1,
                     insert.na.rows = TRUE,
                     subject.gap = NULL,
                     keep.intermediate = FALSE,
                     conduct = TRUE) {
  
  subgroup.time.mode <- match.arg(subgroup.time.mode)
  
  # ---- Argument validation, added in ver 0.0.4 and retained in master ----
  if (is.null(varnames)) stop("varnames must be supplied.")
  if (is.null(dataframe)) stop("dataframe must be supplied.")
  if (is.null(id)) stop("id must be supplied.")
  if (is.null(time)) stop("time must be supplied.")
  if (is.null(directory)) stop("directory must be supplied.")
  if (!is.data.frame(dataframe)) stop("dataframe must be a data.frame.")
  if (!id %in% names(dataframe)) stop("id column not found in dataframe.")
  if (!time %in% names(dataframe)) stop("time column not found in dataframe.")
  missing.vars <- setdiff(varnames, names(dataframe))
  if (length(missing.vars)) {
    stop("Variables not found in dataframe: ", paste(missing.vars, collapse = ", "))
  }
  
  suppressPackageStartupMessages({
    library(parallel)
    library(dynr)
    library(OpenMx)
    library(igraph)
    library(qgraph)
    library(gtools)
  })
  if (isTRUE(conduct) && sub.sig.thrsh != 1.00 && !requireNamespace("nloptr", quietly = TRUE)) {
    stop("Package 'nloptr' is required when conduct = TRUE and subgrouping is enabled.")
  }
  
  # ---- Directory setup ----
  dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(directory, "MIs"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(directory, "Models"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(directory, "Models", "Individuals"), showWarnings = FALSE, recursive = TRUE)
  
  # The OpenMx SSCT model uses a data column named Time.  To preserve the
  # public `time` argument and maintain compatibility with the 0.0.4 script,
  # copy the requested time column into Time before all fitting steps.
  dataframe$Time <- dataframe[[time]]
  time_col <- "Time"
  ids <- unique(dataframe[[id]])
  nvar <- length(varnames)
  if (is.null(PE.var)) PE.var <- diag(1, nvar)
  if (length(ME.var) == 1L) ME.var <- diag(as.numeric(ME.var), nvar)
  if (length(PE.var) == 1L) PE.var <- diag(as.numeric(PE.var), nvar)
  
  if (length(ids) < cores) {
    cores <- length(ids)
    message("Cores adjusted to your sample-size!")
  }
  cores <- max(1, as.integer(cores))
  
  # ---- Local IO/plot helpers ----
  save_rds <- function(object, file) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(pattern = paste0(".", basename(file), "_"), tmpdir = dirname(file))
    saveRDS(object, tmp)
    if (!file.rename(tmp, file)) {
      file.copy(tmp, file, overwrite = TRUE)
      unlink(tmp)
    }
    invisible(file)
  }
  
  safe_png <- function(filename, expr, width = 800, height = 800) {
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    png(filename = filename, width = width, height = height)
    on.exit(dev.off(), add = TRUE)
    tryCatch(force(expr), error = function(e) message("Plot failed for ", filename, ": ", e$message))
    invisible(filename)
  }
  
  get_cells <- function(x) {
    as.numeric(unlist(regmatches(x, gregexpr("\\d+", x))))[1:2]
  }
  
  # ---- Subject-level fitting helper used by group and subgroup shared search ----
  fit_subjects <- function(data_for_fit, subject_ids, DRIFT, model_dir_out,
                           save_models = TRUE, mi_dir_out = file.path(directory, "MIs")) {
    subject_ids <- as.character(subject_ids)
    dir.create(mi_dir_out, recursive = TRUE, showWarnings = FALSE)
    dir.create(model_dir_out, recursive = TRUE, showWarnings = FALSE)
    fit_one <- function(i) {
      subset_dat <- data_for_fit[as.character(data_for_fit[[id]]) == as.character(i), , drop = FALSE]
      if (!nrow(subset_dat)) return(NULL)
      analysis_result <- tryCatch({
        mxTryHard(build_ou_model(subset_dat, DRIFT, nvar, varnames, ME.var, PE.var,
                                 model_name = "OUMod", lb = -10, ub = 10,
                                 time_col = time_col))
      }, error = function(e) {
        message("Error for subject ", i, ": ", e$message)
        NULL
      })
      if (!is.null(analysis_result)) {
        MIs <- tryCatch(JPmx(analysis_result, matrices = "A"), error = function(e) {
          message("MI error for subject ", i, ": ", e$message)
          NULL
        })
        if (!is.null(MIs)) {
          save_rds(MIs, file.path(mi_dir_out, paste0("MI_", i, ".RDS")))
        }
        if (isTRUE(save_models)) {
          save_rds(analysis_result, file.path(model_dir_out, paste0("Model_", i, ".RDS")))
        }
      }
      invisible(analysis_result)
    }
    if (cores <= 1 || length(subject_ids) <= 1) {
      invisible(lapply(subject_ids, fit_one))
    } else {
      cl <- makeCluster(min(cores, length(subject_ids)), type = "PSOCK")
      on.exit(stopCluster(cl), add = TRUE)
      clusterExport(cl, c("fit_one", "JPmx", "build_ou_model", "save_rds"), envir = environment())
      clusterEvalQ(cl, {
        packs <- c("dynr", "OpenMx", "qgraph")
        invisible(lapply(packs, require, character.only = TRUE))
      })
      invisible(parLapply(cl, subject_ids, fit_one))
    }
  }
  
  param_names <- compute_param_names(nvar, model_name = "OUMod")
  
  read_mi_stack <- function(valid_ids = NULL) {
    rdss <- list.files(file.path(directory, "MIs"), pattern = "\\.RDS$", full.names = TRUE)
    if (!is.null(valid_ids)) {
      rdss <- rdss[gsub("MI_|\\.RDS", "", basename(rdss)) %in% as.character(valid_ids)]
    }
    rdss <- mixedsort(rdss)
    agg <- aggregate_mis(rdss, param_names)
    list(files = agg$files, EPCs = agg$EPCs, rdss = rdss)
  }
  
  delete_subject_files <- function(subject_ids, model_dir_out) {
    subject_ids <- as.character(subject_ids)
    unlink(file.path(directory, "MIs", paste0("MI_", subject_ids, ".RDS")), force = TRUE)
    unlink(file.path(model_dir_out, paste0("Model_", subject_ids, ".RDS")), force = TRUE)
  }
  
  prune_shared <- function(DRIFT, model_dir_out, threshold, protected_mask = NULL) {
    rdss1 <- list.files(model_dir_out, pattern = "^Model_.*\\.RDS$", full.names = TRUE)
    if (!length(rdss1)) return(list(DRIFT = DRIFT, pruned = FALSE))
    models <- list()
    for (prn in seq_along(rdss1)) {
      temp1 <- tryCatch(readRDS(rdss1[prn]), error = function(e) {
        message("Failed to read ", prn, ": ", e$message)
        NULL
      })
      if (is.null(temp1)) next
      drifts <- tryCatch(subset(summary(temp1)$parameters, matrix == "A"), error = function(e) NULL)
      if (is.null(drifts) || !nrow(drifts)) next
      temp.mat1 <- matrix(NA, nvar, nvar)
      if ("name" %in% names(drifts)) {
        cells <- matrix(as.numeric(unlist(regmatches(drifts$name, gregexpr("\\d+", drifts$name)))),
                        ncol = 2, byrow = TRUE)
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
          temp.mat1[row_indices[ii], col_indices[ii]] <- ifelse(abs(est) / se > qnorm(0.95), TRUE, FALSE)
        }
      }
      models[[length(models) + 1]] <- cbind(temp.mat1)
    }
    if (!length(models)) return(list(DRIFT = DRIFT, pruned = FALSE))
    arr <- simplify2array(models)
    true_count <- apply(arr, c(1, 2), function(x) sum(x == TRUE, na.rm = TRUE)) / length(models)
    diag(true_count) <- 1.00
    # Higher-level paths are protected from pruning at lower levels.
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
    message(paste0("Pruned drift parameter A[", cells[1], ",", cells[2], "]"))
    list(DRIFT = DRIFT, pruned = TRUE)
  }
  
  shared_search <- function(subject_ids, DRIFT, alpha, threshold,
                            model_dir_out, subgroup = NULL,
                            protected_mask = NULL) {
    iterate <- 0
    count <- 1
    bh_thresholds <- make_threshold_sequence(nvar^2, alpha, ben_hoch = ben.hoch)
    
    while (iterate < 1) {
      stack <- read_mi_stack(valid_ids = subject_ids)
      files <- stack$files
      if (is.null(files) || !ncol(files)) break
      
      # v0.0.5 shared-path selection restores original ctgimme behavior:
      #   1. rank candidate paths by mean MI.Full across usable subjects;
      #   2. evaluate only that mean-best candidate;
      #   3. apply the original adapted BH threshold sequence to individual
      #      p-values for that candidate;
      #   4. stop the current shared search when the mean-best candidate no
      #      longer meets the threshold, rather than scanning to lower-ranked
      #      alternatives.
      files <- abs(files)
      usable_subjects <- colSums(is.finite(files)) > 0
      files <- files[, usable_subjects, drop = FALSE]
      n_converge <- ncol(files)
      n_requested <- length(subject_ids)
      saveRDS(files, paste0(directory, "/LOOKIT.rds"))
      if (is.null(files) || !n_converge) break
      
      if (n_converge <= (n_requested / 2)) {
        message("No shared path added: usable MI files for <= half of subjects.")
        iterate <- 1
        next
      }
      
      has_any_mi <- rowSums(is.finite(files)) > 0
      if (!any(has_any_mi)) break
      
      mean_mi <- rep(-Inf, nrow(files))
      names(mean_mi) <- rownames(files)
      mean_mi[has_any_mi] <- matrixStats::rowMedians(files[has_any_mi, , drop = FALSE], na.rm = TRUE)
      
      best_idx <- which.max(mean_mi)
      if (!length(best_idx) || !is.finite(mean_mi[best_idx])) break
      
      selected_param <- names(mean_mi)[best_idx]
      cells <- get_cells(selected_param)
      if (anyNA(cells)) {
        message("Stopping shared search: malformed mean-best parameter name: ", selected_param)
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
      mean_p <- pchisq(mean_mi[best_idx], df = 1, lower.tail = FALSE)
      
      if (is.finite(sig_prop) && sig_prop >= threshold) {
        DRIFT[cells[1], cells[2]] <- paste0("A_", cells[1], ",", cells[2])
        message(paste0(
          "Adding drift parameter A[", cells[1], ",", cells[2], "]",
          " | mean MI = ", round(mean_mi[best_idx], 4),
          # " | mean-MI p = ", signif(mean_p, 4),
          " | significant in ", sig_count, " of ", n_converge,
          " usable subjects (", round(sig_prop, 3), ")",
          " | adapted BH alpha = ", signif(n_alpha, 4)
        ))
        print(round(selected_mi, 2))
        message(paste0("Completed Step ", count))
        delete_subject_files(subject_ids, model_dir_out)
        fit_subjects(dataframe, subject_ids, DRIFT, model_dir_out, save_models = TRUE)
        count <- count + 1
      } else {
        message(paste0(
          "Stopping shared search: mean-best candidate A[", cells[1], ",", cells[2], "]",
          " did not meet threshold | mean MI = ", round(mean_mi[best_idx], 4),
          " | mean-MI p = ", signif(mean_p, 4),
          " | significant in ", sig_count, " of ", n_converge,
          " usable subjects (", round(sig_prop, 3), ")",
          " | required proportion >= ", threshold,
          " | adapted BH alpha = ", signif(n_alpha, 4)
        ))
        message("IT'S PRUNING TIME!")
        prune <- 0
        while (prune < 1) {
          if (sum(DRIFT != "0") == nvar) {
            prune <- 1
            break
          }
          pruned <- prune_shared(
            DRIFT,
            model_dir_out,
            threshold,
            protected_mask = protected_mask
          )
          DRIFT <- pruned$DRIFT
          if (!isTRUE(pruned$pruned)) {
            prune <- 1
          } else {
            delete_subject_files(subject_ids, model_dir_out)
            fit_subjects(dataframe, subject_ids, DRIFT, model_dir_out, save_models = TRUE)
          }
        }
        iterate <- 1
      }
    }
    DRIFT
  }
  
  build_subgroup_data <- function(new.data) {
    subgroup.data <- NULL
    for (sid in unique(new.data[[id]])) {
      temp <- new.data[as.character(new.data[[id]]) == as.character(sid), , drop = FALSE]
      if (isTRUE(insert.na.rows)) {
        temp <- rbind(temp, rep(NA, ncol(temp)))
        names(temp) <- names(new.data)
      }
      subgroup.data <- rbind(subgroup.data, temp)
    }
    rownames(subgroup.data) <- NULL
    if (subgroup.time.mode %in% c("legacy", "schedule")) {
      nrows <- nrow(subgroup.data)
      cycles <- floor((nrows - 1) / length(measurement.schedule))
      subgroup.data[[time_col]] <- as.vector(
        sapply(0:cycles, function(d) d * cycle.interval + measurement.schedule)
      )[1:nrows]
    }
    if (subgroup.time.mode == "observed") {
      # Preserve observed within-subject timing by default. If subject.gap is
      # supplied, concatenate subjects onto a single increasing time axis.
      if (!is.null(subject.gap)) {
        subgroup.data2 <- NULL
        offset <- 0
        for (sid in unique(new.data[[id]])) {
          nr_sid <- sum(as.character(new.data[[id]]) == as.character(sid)) + as.integer(insert.na.rows)
          temp <- subgroup.data[as.character(subgroup.data[[id]]) == as.character(sid) |
                                  is.na(subgroup.data[[id]]), , drop = FALSE]
          temp <- temp[seq_len(nr_sid), , drop = FALSE]
          temp[[time_col]] <- temp[[time_col]] + offset
          mx <- max(temp[[time_col]], na.rm = TRUE)
          if (is.finite(mx)) offset <- mx + subject.gap
          subgroup.data2 <- rbind(subgroup.data2, temp)
        }
        subgroup.data <- subgroup.data2
        rownames(subgroup.data) <- NULL
      }
    }
    subgroup.data
  }
  
  # ---- Initial scaling and group-level shared search ----
  message("Scaling Variables for Analysis")
  for (ID in ids) {
    idx <- which(as.character(dataframe[[id]]) == as.character(ID))
    dataframe[idx, varnames] <- ctsgimme_safe_scale(dataframe[idx, varnames, drop = FALSE])
  }
  
  unlink(list.files(file.path(directory, "MIs"), pattern = "\\.RDS$", full.names = TRUE), force = TRUE)
  unlink(list.files(file.path(directory, "Models"), pattern = "^Model_.*\\.RDS$", full.names = TRUE), force = TRUE)
  
  DRIFT <- diag(paste0("A_", 1:nvar, ",", 1:nvar), nvar)
  fit_subjects(dataframe, ids, DRIFT, file.path(directory, "Models"), save_models = TRUE)
  DRIFT <- shared_search(ids, DRIFT, Galpha, sig.thrsh, file.path(directory, "Models"))
  G.DRIFT <- DRIFT
  saveRDS(G.DRIFT, file.path(directory, "GStruc.RDS"))
  safe_png(file.path(directory, "Group Paths.png"), {
    qgraph(t(G.DRIFT != "0") * 1, layout = "circle", labels = varnames,
           edge.width = 5, diag = TRUE, edge.labels = "GROUP")
  })
  message("Group Search Complete.")
  
  # ---- Subgroup discovery ----
  walktrap_comm <- NULL
  if (sub.sig.thrsh == 1.00) {
    memb <- rep(1, length(ids))
    names(memb) <- as.character(ids)
    message("Subgrouping Disabled for Testing.")
  } else {
    message("Beginning Subgrouping Stage")
    rdss1 <- mixedsort(list.files(file.path(directory, "Models"), pattern = "^Model_.*\\.RDS$", full.names = TRUE))
    rdss2 <- mixedsort(list.files(file.path(directory, "MIs"), pattern = "^MI_.*\\.RDS$", full.names = TRUE))
    mi_ids <- gsub("MI_|\\.RDS", "", basename(rdss2))
    mi_by_id <- setNames(rdss2, mi_ids)
    models <- list()
    feature.weights <- list()
    kept.ids <- character(0)
    
    # common_path_weight   <- 5.0
    # dominant_path_weight <- 2.5
    # weaker_path_weight   <- 1.0
    # shared_zero_weight   <- 0.25
    common_path_weight   <- 2.0
    dominant_path_weight <- 1.0
    weaker_path_weight   <- 0.50
    shared_zero_weight   <- 0.05
    pair.idx <- which(upper.tri(matrix(NA, nvar, nvar)), arr.ind = TRUE)
    
    for (file in seq_along(rdss1)) {
      model_id <- gsub("Model_|\\.RDS", "", basename(rdss1[file]))
      mi_file <- unname(mi_by_id[model_id])
      if (length(mi_file) == 0 || is.na(mi_file)) next
      temp1 <- tryCatch(readRDS(rdss1[file]), error = function(e) NULL)
      temp2 <- tryCatch(readRDS(mi_file), error = function(e) NULL)
      if (is.null(temp1) || is.null(temp2)) next
      drifts <- tryCatch(subset(summary(temp1)$parameters, matrix == "A"), error = function(e) NULL)
      if (is.null(drifts) || !nrow(drifts)) next
      if ("name" %in% names(drifts)) {
        cells <- matrix(as.numeric(unlist(regmatches(drifts$name, gregexpr("\\d+", drifts$name)))),
                        ncol = 2, byrow = TRUE)
      } else {
        cells <- cbind(as.numeric(drifts$row), as.numeric(drifts$col))
      }
      MI_vals <- temp2$MI.Full
      EPC_vals <- temp2$EPC
      if (is.null(MI_vals) || is.null(EPC_vals) || !length(MI_vals)) next
      MI.cells <- matrix(as.numeric(unlist(regmatches(names(MI_vals), gregexpr("\\d+", names(MI_vals))))),
                         ncol = 2, byrow = TRUE)
      MI.cells <- cbind(MI.cells, MI_vals, EPC_vals)
      
      mi.raw.mat <- matrix(NA_real_, nvar, nvar)
      if (nrow(MI.cells) > 0) {
        for (ii in seq_len(nrow(MI.cells))) {
          r <- MI.cells[ii, 1]
          c <- MI.cells[ii, 2]
          mi.raw <- suppressWarnings(abs(as.numeric(MI.cells[ii, 3])))
          if (is.finite(r) && is.finite(c) && r >= 1 && r <= nvar &&
              c >= 1 && c <= nvar && is.finite(mi.raw)) {
            mi.raw.mat[r, c] <- mi.raw
          }
        }
      }
      MI.cells[, 3] <- ifelse(MI.cells[, 3] > qchisq(0.975, 1), MI.cells[, 4], 0)
      temp.mat1 <- temp.mat2 <- matrix(NA, nvar, nvar)
      temp.w1 <- matrix(common_path_weight, nvar, nvar)
      temp.w2 <- matrix(weaker_path_weight, nvar, nvar)
      
      for (ii in seq_len(nrow(cells))) {
        val <- drifts[ii, "Estimate"]
        se <- drifts[ii, "Std.Error"]
        if (is.finite(val) && is.finite(se) && se > 0) {
          temp.mat1[cells[ii, 1], cells[ii, 2]] <- ifelse(abs(val) / se > qnorm(0.95), val, 0)
        }
      }
      for (ii in seq_len(nrow(MI.cells))) {
        temp.mat2[MI.cells[ii, 1], MI.cells[ii, 2]] <- MI.cells[ii, 3]
      }
      if (nrow(pair.idx) > 0) {
        for (pp in seq_len(nrow(pair.idx))) {
          r <- pair.idx[pp, 1]
          c <- pair.idx[pp, 2]
          mi.forward <- mi.raw.mat[r, c]
          mi.reverse <- mi.raw.mat[c, r]
          if (!is.finite(mi.forward)) mi.forward <- 0
          if (!is.finite(mi.reverse)) mi.reverse <- 0
          if (mi.forward > mi.reverse) {
            temp.w2[r, c] <- dominant_path_weight
            temp.w2[c, r] <- weaker_path_weight
          } else if (mi.reverse > mi.forward) {
            temp.w2[c, r] <- dominant_path_weight
            temp.w2[r, c] <- weaker_path_weight
          } else {
            temp.w2[r, c] <- weaker_path_weight
            temp.w2[c, r] <- weaker_path_weight
          }
        }
      }
      models[[length(models) + 1]] <- cbind(temp.mat1, temp.mat2)
      feature.weights[[length(feature.weights) + 1]] <- cbind(temp.w1, temp.w2)
      kept.ids <- c(kept.ids, model_id)
    }
    
    if (length(models) < 2) {
      memb <- rep(1, length(ids))
      names(memb) <- as.character(ids)
      walktrap_comm <- NULL
      warning("Subgrouping stage had fewer than two valid subject models; assigning all subjects to one subgroup.")
    } else {
      adj.mat <- matrix(NA, length(models), length(models))
      for (ii in seq_along(models)) {
        for (jj in seq_along(models)) {
          x <- c(models[[ii]])
          y <- c(models[[jj]])
          sx <- sign(x)
          sy <- sign(y)
          w <- (c(feature.weights[[ii]]) + c(feature.weights[[jj]])) / 2
          same.nonzero <- !is.na(sx) & !is.na(sy) & sx == sy & sx != 0
          same.zero <- !is.na(sx) & !is.na(sy) & sx == 0 & sy == 0
          nonzero.sim <- sum(w[same.nonzero], na.rm = TRUE)
          zero.sim <- shared_zero_weight * sum(same.zero, na.rm = TRUE)
          adj.mat[ii, jj] <- nonzero.sim + zero.sim
        }
      }
      saveRDS(adj.mat, file = file.path(directory, "adjacencymatrix.RDS"))
      
      jpmodmax1 <- function(x, m) {
        m2 <- m
        threshold_value <- quantile(m2[upper.tri(m2, diag = FALSE)], x, na.rm = TRUE)
        m2[m2 <= threshold_value] <- 0
        g <- igraph::graph_from_adjacency_matrix(m2, mode = "undirected", weighted = TRUE, diag = FALSE)
        p <- igraph::cluster_walktrap(g)
        mem <- igraph::membership(p)
        n <- length(mem)
        sk <- numeric(n)
        ask <- numeric(n)
        for (i in seq_len(n)) {
          internal_i <- 0
          external_i <- 0
          for (j in seq_len(n)) {
            if (m2[i, j] != 0) {
              if (mem[i] == mem[j]) internal_i <- internal_i + m2[i, j]
              else external_i <- external_i + m2[i, j]
            }
          }
          sk[i] <- external_i
          ask[i] <- internal_i + external_i
        }
        temp.1 <- cbind(mem = mem, sk = sk, ask = ask)
        cluster_ids <- sort(unique(mem))
        c_vals <- numeric(length(cluster_ids))
        for (k in seq_along(cluster_ids)) {
          cl <- cluster_ids[k]
          cluster_nodes <- which(mem == cl)
          cut_edges <- sum(temp.1[cluster_nodes, "sk"])
          total_edges <- sum(temp.1[cluster_nodes, "ask"])
          cluster_size <- length(cluster_nodes)
          if (total_edges == 0) c_vals[k] <- 1 else c_vals[k] <- cut_edges / total_edges
          if (cluster_size == 1) c_vals[k] <- max(c_vals[k], 1)
        }
        conductance <- 1 - mean(c_vals)
        -1 * conductance
      }
      if (isTRUE(conduct)) {
        res <- nloptr::nloptr(
          x0 = 0.5, m = adj.mat, eval_f = jpmodmax1,
          lb = .01, ub = .99,
          opts = list("algorithm" = "NLOPT_GN_DIRECT_L", "ftol_rel" = 1.0e-8, maxeval = 1000)
        )
        adj.mat[which(adj.mat <= quantile(adj.mat[upper.tri(adj.mat, diag = FALSE)], res$solution))] <- 0
      } else {
        adj.mat <- as.matrix(adj.mat)
        adj.mat <- adj.mat - min(adj.mat, na.rm = TRUE)
      }
      diag(adj.mat) <- 0
      g <- graph_from_adjacency_matrix(adj.mat, mode = "undirected", weighted = TRUE, diag = FALSE)
      walktrap_comm <- cluster_walktrap(g)
      safe_png(file.path(directory, "walktrap_community_plot.png"), {
        plot(walktrap_comm, g, vertex.size = 15,
             vertex.label = kept.ids,
             main = "Walktrap Community Detection")
      })
      partial_memb <- membership(walktrap_comm)
      names(partial_memb) <- kept.ids
      memb <- rep(1, length(ids))
      names(memb) <- as.character(ids)
      matched <- match(kept.ids, names(memb))
      keep <- !is.na(matched)
      memb[matched[keep]] <- partial_memb[keep]
      if (any(!keep)) {
        warning("Some subgroup model IDs did not match dataframe IDs; unmatched IDs were left in subgroup 1.")
      }
    }
  }
  thewalktrap_comm = walktrap_comm
  print(memb)
  
  # ---- Subgroup-level search and individual fitting ----
  for (subgroup in sort(unique(memb))) {
    DRIFT <- G.DRIFT
    saveRDS(DRIFT, file.path(directory, "GStruc.RDS"))
    subgroup_dir <- file.path(directory, "Models", paste0("Subgroup ", subgroup))
    dir.create(subgroup_dir, showWarnings = FALSE, recursive = TRUE)
    sg.ids <- names(memb)[memb == subgroup]
    new.data <- dataframe[as.character(dataframe[[id]]) %in% as.character(sg.ids), , drop = FALSE]
    valid_ids <- unique(new.data[[id]])
    if (!length(valid_ids)) next
    
    # Protect group-level paths during subgroup-level pruning.
    DRIFT <- shared_search(valid_ids, DRIFT, S.Galpha, sub.sig.thrsh,
                           subgroup_dir, subgroup = subgroup,
                           protected_mask = G.DRIFT != "0")
    message(paste0("Subgroup Search ", subgroup, " of ", max(memb), " Complete."))
    
    if (subgroup.model == TRUE) {
      subgroup.data <- build_subgroup_data(new.data)
      message(paste0("Fitting Parameterized Model of Subgroup ", subgroup))
      fit <- tryCatch(mxTryHard(build_ou_model(subgroup.data, DRIFT, nvar, varnames,
                                               ME.var, PE.var, model_name = "OUMod",
                                               lb = -10, ub = 10, time_col = time_col)),
                      error = function(e) {
                        message("Subgroup parameterized model failed for subgroup ", subgroup, ": ", e$message)
                        NULL
                      })
      if (!is.null(fit)) {
        sum.fit <- summary(fit)
        params.A <- subset(sum.fit$parameters, matrix == "A")
        effects <- matrix(0, nvar, nvar)
        for (ii in seq_len(nrow(params.A))) {
          effects[params.A$col[ii], params.A$row[ii]] <- params.A$Estimate[ii]
        }
        sig <- ifelse(abs(params.A$Estimate) > qnorm(0.95) * params.A$Std.Error, "*", "ns")
        vals <- cbind(round(params.A$Estimate, 2), sig)
        edge_labs <- rep("", nvar * nvar)
        for (ii in seq_len(nrow(params.A))) {
          idx <- (params.A$col[ii] - 1) * nvar + params.A$row[ii]
          edge_labs[idx] <- paste0(vals[ii, 2], " (", vals[ii, 1], ")")
        }
        shared <- G.DRIFT != "0" & DRIFT != "0"
        group_only <- G.DRIFT == "0" & DRIFT != "0"
        colors <- character(length = length(DRIFT))
        colors[shared] <- "gray"
        colors[group_only & effects > 0] <- "blue"
        colors[group_only & effects < 0] <- "red"
        colors[group_only & effects == 0] <- "black"
        colors[DRIFT == "0"] <- "black"
        safe_png(file.path(subgroup_dir, paste0("Subgroup ", subgroup, " Params.png")), {
          qgraph(effects, layout = "circle", labels = varnames,
                 edge.width = 1, diag = TRUE, edge.labels = edge_labs,
                 theme = "colorblind", edge.color = c(colors), fade = TRUE)
        })
        if (!requireNamespace("expm", quietly = TRUE)) {
          stop("Package 'expm' is required for subgroup delta-t plots when subgroup.model = TRUE.")
        }
        for (ints in time.intervals) {
          delt <- round(expm::expm(effects * ints), 3)
          safe_png(file.path(subgroup_dir, paste0("Subgroup ", subgroup, " Delta_t = ", ints, ".png")), {
            qgraph(delt, layout = "circle", labels = varnames, fade = TRUE,
                   edge.width = 1, diag = TRUE, edge.labels = delt, maximum = 1.00,
                   theme = "colorblind", title = paste0("Subgroup ", subgroup, "; Delta_t = ", ints))
          })
        }
        save_rds(fit, file.path(subgroup_dir, paste0("Subgroup_", subgroup, "Model.RDS")))
      }
    }
    
    message("Beginning Individual Model Fitting for Subgroup Members")
    m <- (nvar^2) - sum(DRIFT != "0")
    # Restore original adapted-BH construction for individual additions.
    # Thresholds corresponding to already-free higher-level paths are removed,
    # so individual pruning/addition only operates within the individual level.
    nks <- make_threshold_sequence(nvar^2, Ialpha, ben_hoch = ben.hoch)
    n_higher_paths <- sum(DRIFT != "0")
    if (nrow(nks) > 0 && n_higher_paths > 0) {
      if (n_higher_paths >= nrow(nks)) {
        nks <- matrix(numeric(0), 0, 1)
      } else {
        nks <- nks[-seq_len(n_higher_paths), , drop = FALSE]
      }
    }
    SG.DRIFT <- DRIFT
    saveRDS(SG.DRIFT, file.path(directory, paste0("SGStruc", subgroup, ".RDS")))
    safe_png(file.path(subgroup_dir, paste0("Subgroup ", subgroup, " Paths.png")), {
      qgraph(t(abs(((G.DRIFT != "0") * 1) - ((DRIFT != "0") * 1))),
             layout = "circle", labels = varnames,
             edge.width = 5, diag = TRUE, edge.labels = paste0("SG-", subgroup))
    })
    
    fit_individual <- function(i) {
      DRIFT <- SG.DRIFT
      subset_dat <- new.data[as.character(new.data[[id]]) == as.character(i), , drop = FALSE]
      osc <- build_ou_model(subset_dat, DRIFT, nvar, varnames, ME.var, PE.var,
                            model_name = "OUMod", lb = -10, ub = 10, time_col = time_col)
      fit <- tryCatch(mxTryHard(osc), error = function(e) {
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
          fit <- tryCatch(mxTryHard(osc), error = function(e) NULL)
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
        MIs <- tryCatch(JPmx(fit, matrices = "A"), error = function(e) NULL)
        if (is.null(MIs) || is.null(MIs$MI.Full)) {
          optimization <- 1
          fit <- fit2
        } else if (length(MIs$MI.Full) == length(names(MIs$MI.Full))) {
          max_idx <- which.max(abs(MIs$MI.Full))
          max_val <- abs(MIs$MI.Full[max_idx])
          if (is.finite(max_val) && max_val >= qchisq(1 - nks[count, ], df = 1)) {
            max_name <- names(MIs$MI.Full)[max_idx]
            cells <- get_cells(max_name)
            osc$A$free[cells[1], cells[2]] <- TRUE
            osc$A$labels[cells[1], cells[2]] <- NA_character_
            message(paste0("Adding drift parameter A[", cells[1], ",", cells[2], "]"))
            count <- count + 1
            if (sum(osc$A$free) == nvar^2) optimization <- 1
          } else {
            optimization <- 1
          }
        } else {
          message(paste0("Malformed MI.Full — names and values mismatch. Skipping subject ", i))
          optimization <- 1
          fit <- fit2
        }
      }
      
      message("Pruning Stage.")
      prune <- 0
      prune_iter <- 0
      max_prune_iter <- nvar^2
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
        # Protect both group-level and subgroup-level paths at the individual level.
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
          new_fit <- tryCatch(mxTryHard(osc), error = function(e) NULL)
          if (is.null(new_fit)) {
            osc$A$free[cells[1, 1], cells[1, 2]] <- TRUE
            message("Pruning fit failed; retaining last converged fit.")
            break
          }
          fit <- new_fit
          message(paste0("NOTE: Pruned drift parameter: A[", cells[1, 1], ",", cells[1, 2], "]!"))
        } else {
          message("No Pruning Conducted.")
          prune <- 1
        }
      }
      
      if (!is.null(fit)) {
        save_rds(fit, file.path(directory, "Models", "Individuals", paste0("FinalModel_", i, ".RDS")))
        stat.sig <- subset(summary(fit)$parameters, matrix == "A")
        ests <- matrix(0, nvar, nvar)
        for (jj in seq_len(nrow(stat.sig))) {
          ests[stat.sig$row[jj], stat.sig$col[jj]] <- stat.sig$Estimate[jj]
        }
        ests <- t(ests)
        safe_png(file.path(directory, "Models", "Individuals", paste0("FinalModel_", i, ".PNG")), {
          qgraph(ests, layout = "circle", labels = varnames,
                 edge.width = 1, diag = TRUE, edge.labels = round(ests, 2),
                 theme = "colorblind", fade = FALSE)
        })
      }
      invisible(fit)
    }
    
    if (cores <= 1 || length(valid_ids) <= 1) {
      invisible(lapply(valid_ids, fit_individual))
    } else {
      cl <- makeCluster(min(cores, length(valid_ids)), type = "PSOCK")
      # on.exit(stopCluster(cl), add = TRUE)
      clusterExport(cl, c("fit_individual", "JPmx", "build_ou_model", "save_rds", "safe_png", "get_cells"),
                    envir = environment())
      clusterEvalQ(cl, {
        packs <- c("dynr", "OpenMx", "qgraph")
        invisible(lapply(packs, require, character.only = TRUE))
      })
      invisible(parLapply(cl, valid_ids, fit_individual))
    }
  }
  
  if (!keep.intermediate) {
    unlink(file.path(directory, "MIs"), recursive = TRUE, force = TRUE)
    unlink(list.files(file.path(directory, "Models"), pattern = "^Model_.*\\.RDS$", full.names = TRUE), force = TRUE)
    subgroup_dirs <- list.dirs(file.path(directory, "Models"), recursive = FALSE, full.names = TRUE)
    subgroup_dirs <- subgroup_dirs[grepl("Subgroup ", basename(subgroup_dirs), fixed = TRUE)]
    for (sd in subgroup_dirs) {
      unlink(list.files(sd, pattern = "^Model_.*\\.RDS$", full.names = TRUE), force = TRUE)
    }
  }
  
  message(paste0("Subgrouping with Continuous-Time GIMME Complete. Find networks in ", directory, "."))
  
  if (sub.sig.thrsh == 1) {
    invisible(list(
      message = "Continuous-Time S-GIMME Complete.",
      G.DRIFT = G.DRIFT,
      membership = memb,
      directory = directory
    ))
  } else {
    thewalktrap_comm
  }
}
