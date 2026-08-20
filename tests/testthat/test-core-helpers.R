test_that("within-person scaling is stable for constant columns", {
  input <- cbind(variable = c(1, 2, 3), constant = c(4, 4, 4))
  scaled <- ctgimme:::ctgimme_safe_scale(input)

  expect_equal(as.numeric(scaled[, "variable"]), c(-1, 0, 1))
  expect_equal(as.numeric(scaled[, "constant"]), c(0, 0, 0))
  expect_identical(colnames(scaled), colnames(input))
})

test_that("the OpenMx modification-index API is available", {
  expect_true(is.function(OpenMx::mxMI))
})

test_that("PSOCK worker warnings are relayed without changing return values", {
  cluster <- parallel::makeCluster(2L)
  on.exit(parallel::stopCluster(cluster), add = TRUE)

  worker_function <- function(subject_id) {
    if (identical(subject_id, "B")) warning("synthetic numerical warning")
    if (identical(subject_id, "C")) {
      signalCondition(simpleWarning("directly signaled warning condition"))
    }
    paste0("fit-", subject_id)
  }
  worker_results <- parallel::parLapply(
    cluster,
    c("A", "B", "C"),
    ctgimme:::.ctgimme_capture_worker_result,
    fit_function = worker_function
  )

  warning_messages <- character()
  values <- withCallingHandlers(
    ctgimme:::.ctgimme_relay_worker_results(worker_results),
    warning = function(condition) {
      warning_messages <<- c(warning_messages, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(values, list("fit-A", "fit-B", "fit-C"))
  expect_length(warning_messages, 2L)
  expect_true(any(grepl("subject B", warning_messages, fixed = TRUE)))
  expect_true(any(grepl("synthetic numerical warning", warning_messages,
    fixed = TRUE
  )))
  expect_true(any(grepl("subject C", warning_messages, fixed = TRUE)))
  expect_true(any(grepl("directly signaled warning condition", warning_messages,
    fixed = TRUE
  )))
  expect_silent(
    suppressed <- suppressWarnings(
      ctgimme:::.ctgimme_relay_worker_results(worker_results)
    )
  )
  expect_identical(suppressed, values)

  fatal_function <- function(subject_id) {
    if (identical(subject_id, "B")) stop("synthetic uncaught worker error")
    subject_id
  }
  expect_error(
    parallel::parLapply(
      cluster,
      c("A", "B"),
      ctgimme:::.ctgimme_capture_worker_result,
      fit_function = fatal_function
    ),
    "synthetic uncaught worker error",
    fixed = TRUE
  )
})

test_that("parallel fitting batches relay warnings with subject IDs", {
  output_directory <- tempfile("ctgimme-worker-relay-")
  on.exit(unlink(output_directory, recursive = TRUE), add = TRUE)
  input <- data.frame(
    id = rep(c("A", "B"), each = 2L),
    Time = rep(0:1, 2L),
    x = c(0, 1, 1, 0)
  )
  context <- list(
    id = "id",
    cores = 2L,
    directory = output_directory,
    nvar = 1L,
    varnames = "x",
    ME.var = matrix(0.05, 1L, 1L),
    PE.var = matrix(1, 1L, 1L),
    ME.free = matrix(FALSE, 1L, 1L),
    PE.free = matrix(FALSE, 1L, 1L),
    time_col = "Time",
    verbose = FALSE
  )
  fake_cluster <- structure(list(), class = "ctgimme_test_cluster")

  testthat::local_mocked_bindings(
    build_ou_model = function(...) list(model = TRUE),
    .ctgimme_mx_try_hard = function(...) {
      warning("synthetic shared-fit warning")
      NULL
    },
    .ctgimme_fit_individual = function(context, i, ...) {
      warning("synthetic individual-fit warning")
      paste0("individual-", i)
    },
    clusterExport = function(...) invisible(NULL),
    parLapply = function(cl, values, fun, ...) lapply(values, fun, ...),
    .package = "ctgimme"
  )

  shared_warnings <- character()
  shared_values <- withCallingHandlers(
    ctgimme:::.ctgimme_fit_subjects(
      context,
      input,
      c("A", "B"),
      matrix("A_1,1", 1L, 1L),
      file.path(output_directory, "Models"),
      worker_cluster = fake_cluster
    ),
    warning = function(condition) {
      shared_warnings <<- c(shared_warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(shared_values, list(NULL, NULL))
  expect_length(shared_warnings, 2L)
  expect_true(all(grepl("synthetic shared-fit warning", shared_warnings)))
  expect_true(any(grepl("subject A", shared_warnings, fixed = TRUE)))
  expect_true(any(grepl("subject B", shared_warnings, fixed = TRUE)))

  individual_warnings <- character()
  individual_values <- withCallingHandlers(
    ctgimme:::.ctgimme_run_individual_fits(
      context,
      c("A", "B"),
      input,
      matrix("A_1,1", 1L, 1L),
      0L,
      matrix(numeric(), 0L, 1L),
      0.01,
      worker_cluster = fake_cluster
    ),
    warning = function(condition) {
      individual_warnings <<- c(
        individual_warnings,
        conditionMessage(condition)
      )
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(individual_values, list("individual-A", "individual-B"))
  expect_length(individual_warnings, 2L)
  expect_true(all(grepl("synthetic individual-fit warning", individual_warnings)))
  expect_true(any(grepl("subject A", individual_warnings, fixed = TRUE)))
  expect_true(any(grepl("subject B", individual_warnings, fixed = TRUE)))
})

test_that("parameter names use OpenMx column-major ordering", {
  expect_identical(
    ctgimme:::compute_param_names(2),
    c("OUMod.A[1,1]", "OUMod.A[2,1]", "OUMod.A[1,2]", "OUMod.A[2,2]")
  )
})

test_that("shared-path thresholds retain the adapted BH construction", {
  expect_equal(
    as.numeric(ctgimme:::make_threshold_sequence(4, 0.05)),
    c(0.05, 0.0375, 0.025, 0.0125)
  )
  expect_equal(
    as.numeric(ctgimme:::make_threshold_sequence(3, 0.05, FALSE)),
    rep(0.05, 3)
  )
  expect_equal(nrow(ctgimme:::make_threshold_sequence(0, 0.05)), 0L)
})

test_that("the initial drift has labeled free diagonals and fixed off-diagonals", {
  drift <- ctgimme:::.ctgimme_initial_drift(3)

  expect_identical(diag(drift), c("A_1,1", "A_2,2", "A_3,3"))
  expect_true(all(drift[row(drift) != col(drift)] == "0"))
  expect_false(anyNA(drift))
})

test_that("noise free-parameter specifications are diagonal logical masks", {
  expect_identical(
    ctgimme:::.ctgimme_diagonal_free_mask(TRUE, 3, "PE.free"),
    diag(TRUE, 3)
  )
  expect_identical(
    ctgimme:::.ctgimme_diagonal_free_mask(
      c(TRUE, FALSE, TRUE),
      3,
      "ME.free"
    ),
    diag(c(TRUE, FALSE, TRUE), 3)
  )
  expect_error(
    ctgimme:::.ctgimme_diagonal_free_mask(
      matrix(TRUE, 2, 2),
      2,
      "PE.free"
    ),
    "must not free off-diagonal",
    fixed = TRUE
  )
  expect_error(
    ctgimme:::.ctgimme_diagonal_free_mask(c(1, 0), 2, "ME.free"),
    "must be one logical value",
    fixed = TRUE
  )
})

test_that("noise variance values remain diagonal", {
  expect_identical(
    ctgimme:::.ctgimme_diagonal_variance_values(2, 2, "PE.var"),
    diag(2, 2)
  )
  expect_error(
    ctgimme:::.ctgimme_diagonal_variance_values(
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
    id = rep("S1", 6),
    Time = 0:5,
    x = c(1, 2, 4, 3, 5, 7),
    y = c(2, 1, 3, 5, 4, 8)
  )
  drift <- matrix(c("A_1,1", "0", "0", "A_2,2"), 2, 2)
  model <- ctgimme:::build_ou_model(
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
  expect_identical(names(model@data@observed), c("Time", "x", "y"))
})

test_that("measurement and process noise diagonals can be selectively free", {
  input <- data.frame(
    Time = 0:5,
    x = c(1, 2, 4, 3, 5, 7),
    y = c(2, 1, 3, 5, 4, 8)
  )
  drift <- matrix(c("A_1,1", "0", "0", "A_2,2"), 2, 2)
  model <- ctgimme:::build_ou_model(
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

test_that("every internal OpenMx retry fit honors the verbose setting", {
  count_fixed <- function(text, pattern) {
    matches <- gregexpr(pattern, text, fixed = TRUE)[[1L]]
    if (length(matches) == 1L && matches[[1L]] == -1L) 0L else length(matches)
  }
  fitters <- list(
    ctgimme:::.ctgimme_fit_subjects,
    ctgimme:::.ctgimme_fit_individual,
    ctgimme:::.ctgimme_fit_subgroup_model
  )

  for (fitter in fitters) {
    code <- paste(deparse(body(fitter)), collapse = "\n")
    expect_gt(count_fixed(code, ".ctgimme_mx_try_hard"), 0L)
    expect_equal(count_fixed(code, "OpenMx::mxTryHard"), 0L)
  }

  helper_code <- paste(
    deparse(body(ctgimme:::.ctgimme_mx_try_hard)),
    collapse = "\n"
  )
  expect_match(helper_code, "optimizer_silent <- FALSE", fixed = TRUE)
  expect_match(helper_code, "suppressMessages", fixed = TRUE)
  expect_match(helper_code, "utils::capture.output", fixed = TRUE)
})

test_that("subgroup edge labels are square matrices", {
  labels <- paste0("path-", seq_len(36))
  label_matrix <- ctgimme:::.ctgimme_edge_label_matrix(labels, 6)

  expect_identical(dim(label_matrix), c(6L, 6L))
  expect_identical(as.vector(label_matrix), labels)
  expect_error(
    ctgimme:::.ctgimme_edge_label_matrix(labels[-1], 6),
    "exactly nvar^2",
    fixed = TRUE
  )
})

test_that("failed plots do not leave blank PNG artifacts", {
  filename <- tempfile(fileext = ".png")

  expect_warning(
    result <- ctgimme:::.ctgimme_safe_png(
      filename,
      stop("intentional plot failure")
    ),
    "intentional plot failure",
    fixed = TRUE
  )
  expect_null(result)
  expect_false(file.exists(filename))
})

test_that("artifact subject IDs preserve identifier text", {
  expect_identical(
    ctgimme:::.ctgimme_artifact_subject_id(
      c("MI_MI_A.RDS.RDS", "MI_Model_B.RDS", "not-an-artifact"),
      "MI_"
    ),
    c("MI_A.RDS", "Model_B", NA_character_)
  )
  expect_identical(
    ctgimme:::.ctgimme_artifact_subject_id(
      "Model_Model_A.RDS.RDS",
      "Model_"
    ),
    "Model_A.RDS"
  )
})
