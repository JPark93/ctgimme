test_that("within-person scaling is stable for constant columns", {
  input <- cbind(variable = c(1, 2, 3), constant = c(4, 4, 4))
  scaled <- ctsgimme:::ctsgimme_safe_scale(input)

  expect_equal(as.numeric(scaled[, "variable"]), c(-1, 0, 1))
  expect_equal(as.numeric(scaled[, "constant"]), c(0, 0, 0))
  expect_identical(colnames(scaled), colnames(input))
})

test_that("the OpenMx modification-index API is available", {
  expect_true(is.function(OpenMx::mxMI))
})

test_that("parameter names use OpenMx column-major ordering", {
  expect_identical(
    ctsgimme:::compute_param_names(2),
    c("OUMod.A[1,1]", "OUMod.A[2,1]", "OUMod.A[1,2]", "OUMod.A[2,2]")
  )
})

test_that("shared-path thresholds retain the adapted BH construction", {
  expect_equal(
    as.numeric(ctsgimme:::make_threshold_sequence(4, 0.05)),
    c(0.05, 0.0375, 0.025, 0.0125)
  )
  expect_equal(
    as.numeric(ctsgimme:::make_threshold_sequence(3, 0.05, FALSE)),
    rep(0.05, 3)
  )
  expect_equal(nrow(ctsgimme:::make_threshold_sequence(0, 0.05)), 0L)
})

test_that("the initial drift has labeled free diagonals and fixed off-diagonals", {
  drift <- ctsgimme:::.ctsgimme_initial_drift(3)

  expect_identical(diag(drift), c("A_1,1", "A_2,2", "A_3,3"))
  expect_true(all(drift[row(drift) != col(drift)] == "0"))
  expect_false(anyNA(drift))
})

test_that("noise free-parameter specifications are diagonal logical masks", {
  expect_identical(
    ctsgimme:::.ctsgimme_diagonal_free_mask(TRUE, 3, "PE.free"),
    diag(TRUE, 3)
  )
  expect_identical(
    ctsgimme:::.ctsgimme_diagonal_free_mask(
      c(TRUE, FALSE, TRUE),
      3,
      "ME.free"
    ),
    diag(c(TRUE, FALSE, TRUE), 3)
  )
  expect_error(
    ctsgimme:::.ctsgimme_diagonal_free_mask(
      matrix(TRUE, 2, 2),
      2,
      "PE.free"
    ),
    "must not free off-diagonal",
    fixed = TRUE
  )
  expect_error(
    ctsgimme:::.ctsgimme_diagonal_free_mask(c(1, 0), 2, "ME.free"),
    "must be one logical value",
    fixed = TRUE
  )
})

test_that("noise variance values remain diagonal", {
  expect_identical(
    ctsgimme:::.ctsgimme_diagonal_variance_values(2, 2, "PE.var"),
    diag(2, 2)
  )
  expect_error(
    ctsgimme:::.ctsgimme_diagonal_variance_values(
      matrix(c(1, 0.2, 0.2, 1), 2, 2),
      2,
      "PE.var"
    ),
    "full covariances are not supported",
    fixed = TRUE
  )
})

test_that("the empirical initial covariance is symmetric and positive definite", {
  input <- data.frame(
    Time = 0:5,
    x = c(1, 2, 4, 3, 5, 7),
    y = c(2, 1, 3, 5, 4, 8)
  )
  drift <- matrix(c("A_1,1", "0", "0", "A_2,2"), 2, 2)
  model <- ctsgimme:::build_ou_model(
    input,
    drift,
    nvar = 2,
    varnames = c("x", "y"),
    ME.var = diag(1e-5, 2),
    PE.var = diag(1, 2)
  )
  p0 <- as.matrix(model$P0$values)

  expect_equal(p0, t(p0), tolerance = 1e-12)
  expect_true(all(is.finite(p0)))
  expect_gt(min(eigen(p0, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_true(any(abs(p0[row(p0) != col(p0)]) > 0))
})

test_that("measurement and process noise diagonals can be selectively free", {
  input <- data.frame(
    Time = 0:5,
    x = c(1, 2, 4, 3, 5, 7),
    y = c(2, 1, 3, 5, 4, 8)
  )
  drift <- matrix(c("A_1,1", "0", "0", "A_2,2"), 2, 2)
  model <- ctsgimme:::build_ou_model(
    input,
    drift,
    nvar = 2,
    varnames = c("x", "y"),
    ME.var = diag(c(0.1, 0.2), 2),
    PE.var = diag(c(0.8, 1.2), 2),
    ME.free = diag(c(TRUE, FALSE), 2),
    PE.free = diag(TRUE, 2)
  )

  expect_identical(diag(model$R$free), c(TRUE, FALSE))
  expect_identical(diag(model$Q$free), c(TRUE, TRUE))
  expect_false(any(model$R$free[row(model$R$free) != col(model$R$free)]))
  expect_false(any(model$Q$free[row(model$Q$free) != col(model$Q$free)]))
  expect_equal(diag(model$R$values), c(0.1, 0.2))
  expect_equal(diag(model$Q$values), c(0.8, 1.2))
})

test_that("every internal OpenMx retry fit avoids the interactive callback", {
  count_fixed <- function(text, pattern) {
    matches <- gregexpr(pattern, text, fixed = TRUE)[[1L]]
    if (identical(matches, -1L)) 0L else length(matches)
  }
  fitters <- list(
    ctsgimme:::fit_and_save_worker,
    ctsgimme:::.ctsgimme_fit_subjects,
    ctsgimme:::.ctsgimme_fit_individual,
    ctsgimme:::.ctsgimme_fit_subgroup_model
  )

  for (fitter in fitters) {
    code <- paste(deparse(body(fitter)), collapse = "\n")
    expect_gt(count_fixed(code, "mxTryHard"), 0L)
    expect_equal(
      count_fixed(code, "mxTryHard"),
      count_fixed(code, "silent = FALSE")
    )
  }
})

test_that("subgroup edge labels are square matrices", {
  labels <- paste0("path-", seq_len(36))
  label_matrix <- ctsgimme:::.ctsgimme_edge_label_matrix(labels, 6)

  expect_identical(dim(label_matrix), c(6L, 6L))
  expect_identical(as.vector(label_matrix), labels)
  expect_error(
    ctsgimme:::.ctsgimme_edge_label_matrix(labels[-1], 6),
    "exactly nvar^2",
    fixed = TRUE
  )
})

test_that("failed plots do not leave blank PNG artifacts", {
  filename <- tempfile(fileext = ".png")

  expect_message(
    result <- ctsgimme:::.ctsgimme_safe_png(
      filename,
      stop("intentional plot failure")
    ),
    "intentional plot failure",
    fixed = TRUE
  )
  expect_null(result)
  expect_false(file.exists(filename))
})
