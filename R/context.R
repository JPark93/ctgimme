# ============================================================================
# Input validation, run context, and file helpers
# ============================================================================

.ctsgimme_validate_inputs <- function(varnames, dataframe, id, time, directory) {
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
  invisible(TRUE)
}

.ctsgimme_validate_subgroup_options <- function(
    sub.sig.thrsh, subgroup.method, max.subgroups) {
  if (sub.sig.thrsh == 1.00 || !identical(subgroup.method, "pam")) {
    return(invisible(TRUE))
  }
  valid.maximum <- length(max.subgroups) == 1L &&
    is.numeric(max.subgroups) &&
    !is.na(max.subgroups) &&
    is.finite(max.subgroups) &&
    max.subgroups >= 2 &&
    max.subgroups <= .Machine$integer.max &&
    max.subgroups == floor(max.subgroups)
  if (!valid.maximum) {
    stop(
      "max.subgroups must be supplied as one integer >= 2 when ",
      "subgroup.method = 'pam'."
    )
  }
  invisible(TRUE)
}

.ctsgimme_validate_time_intervals <- function(time.intervals) {
  valid <- is.numeric(time.intervals) &&
    length(time.intervals) >= 1L &&
    !anyNA(time.intervals) &&
    all(is.finite(time.intervals)) &&
    all(time.intervals >= 0)
  if (!valid) {
    stop("time.intervals must be a nonempty numeric vector of finite, nonnegative values.")
  }
  invisible(TRUE)
}

.ctsgimme_resolve_cores <- function(cores, subject_count) {
  valid_cores <- length(cores) == 1L &&
    is.numeric(cores) &&
    !is.na(cores) &&
    is.finite(cores) &&
    cores >= 1 &&
    cores <= .Machine$integer.max &&
    cores == floor(cores)
  if (!valid_cores) {
    stop("cores must be supplied as one positive integer.")
  }

  requested <- as.integer(cores)
  resolved <- min(requested, as.integer(subject_count))
  if (resolved < requested) {
    message("Cores adjusted to the number of subjects.")
  }

  check_limit <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = ""))
  limit_during_check <- nzchar(check_limit) && check_limit != "false"
  if (limit_during_check && resolved > 2L) {
    resolved <- 2L
    message("Cores adjusted to two for CRAN-style package checking.")
  }

  max(1L, resolved)
}

.ctsgimme_worker_registry <- local({
  registry <- new.env(parent = emptyenv())
  registry$next_id <- 0L
  registry
})

.ctsgimme_register_worker_cluster <- function(cluster) {
  pids <- as.integer(unlist(parallel::clusterCall(cluster, Sys.getpid)))
  handles <- lapply(pids, ps::ps_handle)

  .ctsgimme_worker_registry$next_id <-
    .ctsgimme_worker_registry$next_id + 1L
  key <- paste(
    Sys.getpid(),
    .ctsgimme_worker_registry$next_id,
    sep = ":"
  )
  assign(
    key,
    list(pids = pids, handles = handles),
    envir = .ctsgimme_worker_registry
  )
  attr(cluster, "ctsgimme.worker.key") <- key
  cluster
}

.ctsgimme_make_worker_cluster <- function(cores) {
  cores <- as.integer(cores)
  if (cores <= 1L) return(NULL)

  # Package loading is one of the largest fixed costs on Windows. Keep these
  # PSOCK sessions alive across group, subgroup, and individual fitting batches
  # so that initialization is paid once per CTSGIMME analysis.
  cluster <- makeCluster(cores, type = "PSOCK")
  initialized <- FALSE
  on.exit({
    if (!initialized) {
      .ctsgimme_stop_worker_cluster(cluster)
    }
  }, add = TRUE)

  cluster <- .ctsgimme_register_worker_cluster(cluster)

  status <- clusterEvalQ(cluster, {
    packages <- c("OpenMx", "qgraph")
    loaded <- vapply(
      packages,
      require,
      logical(1),
      character.only = TRUE,
      quietly = TRUE
    )
    if (!all(loaded)) {
      stop(
        "Failed to load worker packages: ",
        paste(packages[!loaded], collapse = ", ")
      )
    }
    OpenMx::mxOption(NULL, "Number of Threads", 1L)
    TRUE
  })
  if (!all(vapply(status, isTRUE, logical(1)))) {
    stop("One or more CTSGIMME workers failed to initialize.")
  }

  initialized <- TRUE
  cluster
}

.ctsgimme_stop_worker_cluster <- function(cluster) {
  if (is.null(cluster)) return(invisible(NULL))

  key <- attr(cluster, "ctsgimme.worker.key", exact = TRUE)
  record <- if (!is.null(key) &&
      exists(key, envir = .ctsgimme_worker_registry, inherits = FALSE)) {
    get(key, envir = .ctsgimme_worker_registry, inherits = FALSE)
  } else {
    NULL
  }

  graceful_error <- NULL
  tryCatch(
    stopCluster(cluster),
    error = function(error) {
      graceful_error <<- error
      invisible(NULL)
    }
  )

  forced <- FALSE
  all_stopped <- TRUE
  if (!is.null(record) && length(record$handles)) {
    exited <- tryCatch(
      ps::ps_wait(record$handles, timeout = 2000L),
      error = function(error) rep(FALSE, length(record$handles))
    )
    survivors <- record$handles[!exited]

    if (length(survivors)) {
      forced <- TRUE
      tryCatch(
        ps::ps_kill(survivors, grace = 200L),
        error = function(error) invisible(NULL)
      )
      stopped_after_kill <- tryCatch(
        ps::ps_wait(survivors, timeout = 5000L),
        error = function(error) rep(FALSE, length(survivors))
      )
      all_stopped <- all(stopped_after_kill)
    }
  }

  if (!is.null(key) &&
      exists(key, envir = .ctsgimme_worker_registry, inherits = FALSE)) {
    rm(list = key, envir = .ctsgimme_worker_registry)
  }

  if (!all_stopped) {
    tryCatch(
      warning(
        "One or more CTSGIMME workers remained alive after forced shutdown.",
        call. = FALSE
      ),
      error = function(error) invisible(NULL)
    )
  } else if (!is.null(graceful_error) || forced) {
    tryCatch(
      warning(
        "Graceful CTSGIMME worker shutdown failed; surviving workers were ",
        "force-terminated.",
        call. = FALSE
      ),
      error = function(error) invisible(NULL)
    )
  }

  invisible(all_stopped)
}

.ctsgimme_load_dependencies <- function(
    conduct, sub.sig.thrsh, subgroup.method = "pam", subgroup.model = FALSE) {
  if (sub.sig.thrsh != 1.00 && identical(subgroup.method, "pam") &&
      !requireNamespace("cluster", quietly = TRUE)) {
    stop("Package 'cluster' is required for PAM subgrouping.")
  }
  if (sub.sig.thrsh != 1.00 && identical(subgroup.method, "legacy") &&
      !requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for legacy subgrouping.")
  }
  if (isTRUE(conduct) && sub.sig.thrsh != 1.00 &&
      identical(subgroup.method, "legacy") &&
      !requireNamespace("nloptr", quietly = TRUE)) {
    stop(
      "Package 'nloptr' is required when conduct = TRUE and ",
      "subgroup.method = 'legacy'."
    )
  }
  if (isTRUE(subgroup.model) && !requireNamespace("expm", quietly = TRUE)) {
    stop("Package 'expm' is required when subgroup.model = TRUE.")
  }
  invisible(TRUE)
}

.ctsgimme_diagonal_free_mask <- function(mask, nvar, argument) {
  valid_scalar <- is.logical(mask) && length(mask) == 1L && !is.na(mask)
  if (valid_scalar) {
    return(diag(mask, nvar))
  }

  valid_vector <- is.logical(mask) && is.null(dim(mask)) &&
    length(mask) == nvar && !anyNA(mask)
  if (valid_vector) {
    return(diag(mask, nvar))
  }

  valid_matrix <- is.logical(mask) && is.matrix(mask) &&
    all(dim(mask) == c(nvar, nvar)) && !anyNA(mask)
  if (valid_matrix) {
    off_diagonal <- row(mask) != col(mask)
    if (any(mask[off_diagonal])) {
      stop(argument, " must not free off-diagonal covariance parameters.")
    }
    return(mask)
  }

  stop(
    argument,
    " must be one logical value, a logical vector with one value per ",
    "process variable, or a diagonal logical matrix."
  )
}

.ctsgimme_diagonal_variance_values <- function(
    values, nvar, argument, default = NULL) {
  if (is.null(values) && !is.null(default)) {
    values <- default
  }
  if (!is.numeric(values) || anyNA(values) || any(!is.finite(values))) {
    stop(argument, " must contain finite numeric variance values.")
  }
  if (length(values) == 1L) {
    values <- diag(as.numeric(values), nvar)
  } else if (!is.matrix(values) ||
             !all(dim(values) == c(nvar, nvar))) {
    stop(
      argument,
      " must be one numeric value or a diagonal numeric matrix with one ",
      "row and column per process variable."
    )
  }

  off_diagonal <- row(values) != col(values)
  if (any(abs(values[off_diagonal]) > sqrt(.Machine$double.eps))) {
    stop(argument, " must be diagonal; full covariances are not supported.")
  }
  diagonal <- diag(values)
  if (any(diagonal < 0)) {
    stop(argument, " must not contain negative variances.")
  }
  diag(diagonal, nvar)
}

.ctsgimme_prepare_context <- function(
    varnames, dataframe, id, time, cores, directory, ME.var, PE.var,
    ben.hoch, ME.free = FALSE, PE.free = FALSE) {
  nvar <- length(varnames)
  ME.var <- .ctsgimme_diagonal_variance_values(ME.var, nvar, "ME.var")
  PE.var <- .ctsgimme_diagonal_variance_values(
    PE.var,
    nvar,
    "PE.var",
    default = diag(1, nvar)
  )
  ME.free <- .ctsgimme_diagonal_free_mask(ME.free, nvar, "ME.free")
  PE.free <- .ctsgimme_diagonal_free_mask(PE.free, nvar, "PE.free")
  if (any(diag(ME.var)[diag(ME.free)] < 1e-8)) {
    stop("Free ME.var starting values must be at least 1e-8.")
  }
  if (any(diag(PE.var)[diag(PE.free)] < 1e-5)) {
    stop("Free PE.var starting values must be at least 1e-5.")
  }

  dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(directory, "MIs"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(directory, "Models"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(directory, "Models", "Individuals"),
             showWarnings = FALSE, recursive = TRUE)

  dataframe$Time <- dataframe[[time]]
  time_col <- "Time"
  ids <- unique(dataframe[[id]])

  cores <- .ctsgimme_resolve_cores(cores, length(ids))

  list(
    varnames = varnames,
    dataframe = dataframe,
    id = id,
    time = time,
    time_col = time_col,
    cores = cores,
    directory = directory,
    ids = ids,
    nvar = nvar,
    ME.var = ME.var,
    PE.var = PE.var,
    ME.free = ME.free,
    PE.free = PE.free,
    ben.hoch = ben.hoch,
    param_names = compute_param_names(nvar, model_name = "OUMod")
  )
}

.ctsgimme_save_rds <- function(object, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(file), "_"), tmpdir = dirname(file))
  saveRDS(object, tmp)
  if (!file.rename(tmp, file)) {
    file.copy(tmp, file, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(file)
}

.ctsgimme_safe_png <- function(filename, expr, width = 800, height = 800) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  device <- NULL
  success <- FALSE
  tryCatch(
    {
      grDevices::png(filename = filename, width = width, height = height)
      device <- grDevices::dev.cur()
      force(expr)
      success <- TRUE
    },
    error = function(e) {
      message("Plot failed for ", filename, ": ", e$message)
    },
    finally = {
      if (!is.null(device) && device %in% grDevices::dev.list()) {
        grDevices::dev.off(device)
      }
      if (!success && file.exists(filename)) {
        unlink(filename)
      }
    }
  )
  invisible(if (success) filename else NULL)
}

.ctsgimme_get_cells <- function(x) {
  as.numeric(unlist(regmatches(x, gregexpr("\\d+", x))))[1:2]
}
