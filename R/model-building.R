# Read a named vector safely from a saved model
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
    warning(
      "Failed to read ", element, " from ", file, ": ", e$message,
      call. = FALSE
    )
    return(NULL)
  })

  if (!is.null(vec)) {
    intersect_names = intersect(names(vec), param_names)
    out[intersect_names] = vec[intersect_names]
  }

  return(out)
}


# Scale subject-level data robustly
ctgimme_safe_scale <- function(x) {
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

# Construct OpenMx parameter names
compute_param_names <- function(nvar, model_name = "OUMod") {
  param_names = character(0)
  for (j in 1:nvar) {
    for (i in 1:nvar) {
      param_names = c(param_names, sprintf("%s.A[%d,%d]", model_name, i, j))
    }
  }
  return(param_names)
}

# Extract matrix indices from parameter names
extract_matrix_indices <- function(names_vec) {
  nums = as.numeric(unlist(regmatches(names_vec, gregexpr("\\d+", names_vec))))
  matrix(nums, ncol = 2, byrow = TRUE)
}

# Construct adapted Benjamini-Hochberg thresholds
make_threshold_sequence <- function(m, alpha, ben_hoch = TRUE) {
  # The search consumes (k / m) * alpha in decreasing order, beginning at
  # alpha and becoming progressively stricter.
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

.ctgimme_initial_covariance <- function(data, varnames, ME.var) {
  observed = as.matrix(data[, varnames, drop = FALSE])
  p0_values = if (nrow(observed) >= 2L) {
    empirical = stats::cov(observed, use = "pairwise.complete.obs")
    latent = empirical - as.matrix(ME.var)
    latent = (latent + t(latent)) / 2
    if (all(is.finite(latent))) {
      # Preserve marginal variances while shrinking noisy cross-covariances.
      latent = 0.5 * latent + 0.5 * diag(diag(latent), nrow(latent))
      eig = eigen(latent, symmetric = TRUE)
      positive_diag = diag(latent)[is.finite(diag(latent)) & diag(latent) > 0]
      scale_floor = if (length(positive_diag)) stats::median(positive_diag) else 1
      eigen_floor = max(1e-8, 1e-6 * scale_floor)
      eig$vectors %*% diag(pmax(eig$values, eigen_floor), nrow(latent)) %*%
        t(eig$vectors)
    } else {
      diag(1, nrow = ncol(observed))
    }
  } else {
    diag(1, nrow = ncol(observed))
  }
  (p0_values + t(p0_values)) / 2
}

.ctgimme_resolve_initial_covariance <- function(
    data, varnames, ME.var, P0.values = NULL) {
  if (is.null(P0.values)) {
    return(.ctgimme_initial_covariance(data, varnames, ME.var))
  }
  P0.values = as.matrix(P0.values)
  nvar = length(varnames)
  if (!identical(dim(P0.values), c(nvar, nvar))) {
    stop("P0.values must have one row and column per process variable.")
  }
  if (any(!is.finite(P0.values))) {
    stop("P0.values must contain only finite values.")
  }
  P0.values = (P0.values + t(P0.values)) / 2
  if (min(eigen(P0.values, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("P0.values must be positive definite.")
  }
  P0.values
}

# Construct a single-series Ornstein-Uhlenbeck model
build_ou_model <- function(data, drift_matrix, nvar, varnames, ME.var, PE.var,
                           ME.free = FALSE, PE.free = FALSE,
                           model_name = "OUMod", lb = -Inf, ub = Inf,
                           time_col = "Time", P0.values = NULL) {
  model_data = data[, c(time_col, varnames), drop = FALSE]
  amat = mxMatrix("Full", nvar, nvar,
                  free   = drift_matrix != "0",
                  name   = "A",
                  ubound = ub, lbound = lb)
  bmat = mxMatrix('Zero', nvar, nvar, name='B')
  cdim = list(varnames, paste0('F', 1:nvar))
  cmat = mxMatrix('Diag', nvar, nvar, FALSE, 1, name = 'C', dimnames = cdim)
  dmat = mxMatrix('Zero', nvar, nvar, name='D')
  qmat = mxMatrix('Diag', nvar, nvar, PE.free, PE.var, name='Q', lbound = 1e-5)
  rmat = mxMatrix('Diag', nvar, nvar, ME.free, ME.var, name='R', lbound = 0)
  xmat = mxMatrix('Full', nvar, 1, FALSE, rep(0, nvar), name='x0')
  p0_values = .ctgimme_resolve_initial_covariance(
    data,
    varnames,
    ME.var,
    P0.values
  )
  pmat = mxMatrix('Full', nvar, nvar, FALSE, p0_values, name='P0')
  umat = mxMatrix('Zero', nvar, 1, name='u')
  tmat = mxMatrix('Full', 1, 1, FALSE, name='time', labels=paste0('data.', time_col))
  mxModel(model_name,
          amat, bmat, cmat, dmat, qmat,
          rmat, xmat, pmat, umat, tmat,
          mxExpectationSSCT('A', 'B', 'C', 'D', 'Q',
                            'R', 'x0', 'P0', 'u', 'time'),
          mxFitFunctionML(),
          mxData(model_data, 'raw'))
}

# Construct one shared-parameter OU model with an independent filter per subject.
build_multisubject_ou_model <- function(
    data, id_col, drift_matrix, nvar, varnames, ME.var, PE.var,
    ME.free = FALSE, PE.free = FALSE, model_name = "OUMod",
    lb = -Inf, ub = Inf, time_col = "Time", P0.values = NULL) {
  if (!id_col %in% names(data)) stop("id_col was not found in data.")
  if (!time_col %in% names(data)) stop("time_col was not found in data.")
  if (anyNA(data[[id_col]])) stop("Subject identifiers must not be missing.")

  ids = unique(data[[id_col]])
  if (!length(ids)) stop("At least one subject is required.")

  free_drift = drift_matrix != "0"
  drift_values = matrix(0, nrow = nvar, ncol = nvar)
  diag(drift_values) = ifelse(diag(free_drift), -0.5, 0)
  amat = mxMatrix(
    "Full", nvar, nvar,
    free = free_drift,
    values = drift_values,
    name = "A",
    ubound = ub,
    lbound = lb
  )
  qmat = mxMatrix(
    "Diag", nvar, nvar, PE.free, PE.var,
    name = "Q", lbound = 1e-5
  )
  rmat = mxMatrix(
    "Diag", nvar, nvar, ME.free, ME.var,
    name = "R", lbound = 0
  )
  bmat = mxMatrix("Zero", nvar, nvar, name = "B")
  cmat = mxMatrix(
    "Diag", nvar, nvar, FALSE, 1,
    name = "C",
    dimnames = list(varnames, paste0("F", seq_len(nvar)))
  )
  dmat = mxMatrix("Zero", nvar, nvar, name = "D")
  xmat = mxMatrix("Full", nvar, 1, FALSE, rep(0, nvar), name = "x0")
  umat = mxMatrix("Zero", nvar, 1, name = "u")

  # Use one common fixed initialization for all independent filters. This
  # preserves the one-subgroup-model interpretation and makes the only
  # subject-specific element the observed trajectory itself.
  p0_values = .ctgimme_resolve_initial_covariance(
    data,
    varnames,
    ME.var,
    P0.values
  )
  pmat = mxMatrix("Full", nvar, nvar, FALSE, p0_values, name = "P0")

  child_names = sprintf("subject_%06d", seq_along(ids))
  children = lapply(seq_along(ids), function(index) {
    subject_data = data[
      as.character(data[[id_col]]) == as.character(ids[[index]]),
      c(varnames, time_col),
      drop = FALSE
    ]
    if (!nrow(subject_data)) stop("A subject likelihood block had no rows.")
    subject_time = subject_data[[time_col]]
    if (!is.numeric(subject_time) || anyNA(subject_time) ||
        any(!is.finite(subject_time))) {
      stop("Each subject's observation times must be finite numeric values.")
    }
    if (is.unsorted(subject_time)) {
      stop("Rows must be ordered by time within each subject.")
    }

    # x0 and P0 initialize each independent filter at its own local time zero.
    # Rebasing removes arbitrary calendar offsets without changing any elapsed
    # time between observations.
    subject_data[[time_col]] = subject_time - subject_time[[1L]]

    mxModel(
      child_names[[index]],
      mxMatrix(
        "Full", 1, 1, FALSE,
        name = "time",
        labels = paste0("data.", time_col)
      ),
      mxExpectationSSCT(
        paste0(model_name, ".A"),
        paste0(model_name, ".B"),
        paste0(model_name, ".C"),
        paste0(model_name, ".D"),
        paste0(model_name, ".Q"),
        paste0(model_name, ".R"),
        paste0(model_name, ".x0"),
        paste0(model_name, ".P0"),
        paste0(model_name, ".u"),
        "time"
      ),
      mxFitFunctionML(),
      mxData(subject_data, "raw")
    )
  })

  model = mxModel(
    model_name,
    amat,
    bmat,
    cmat,
    dmat,
    qmat,
    rmat,
    xmat,
    pmat,
    umat,
    children,
    mxFitFunctionMultigroup(paste0(child_names, ".fitfunction"))
  )
  attr(model, "ctgimme.subject.ids") = as.character(ids)
  attr(model, "ctgimme.subgroup.likelihood") = "multisubject"
  model
}
