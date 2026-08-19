multisubject_test_data <- function(reverse_subjects = FALSE) {
  ids <- if (reverse_subjects) c("B", "A") else c("A", "B")
  values <- list(
    A = c(0.8, 0.5, 0.3, 0.2),
    B = c(-0.7, -0.4, -0.2, -0.1)
  )
  do.call(rbind, lapply(ids, function(id) {
    start <- if (id == "A") 0 else 10
    data.frame(id = id, x = values[[id]], Time = start + c(0, 1, 2, 3))
  }))
}

test_that("multisubject OU uses one shared parameter vector", {
  input <- multisubject_test_data()
  model <- ctsgimme:::build_multisubject_ou_model(
    input,
    id_col = "id",
    drift_matrix = matrix("A_1,1", 1, 1),
    nvar = 1,
    varnames = "x",
    ME.var = matrix(0.1, 1, 1),
    PE.var = matrix(1, 1, 1),
    ME.free = TRUE,
    PE.free = TRUE,
    P0.values = matrix(1, 1, 1)
  )

  expect_s4_class(model, "MxModel")
  expect_length(model@submodels, 2L)
  expect_identical(
    names(OpenMx::omxGetParameters(model)),
    c("OUMod.A[1,1]", "OUMod.Q[1,1]", "OUMod.R[1,1]")
  )
  expect_identical(attr(model, "ctsgimme.subject.ids"), c("A", "B"))
  expect_identical(
    attr(model, "ctsgimme.subgroup.likelihood"),
    "multisubject"
  )
  expect_equal(as.matrix(model$P0$values), matrix(1, 1, 1))

  child_rows <- vapply(
    model@submodels,
    function(child) nrow(child@data@observed),
    integer(1)
  )
  expect_identical(child_rows, c(subject_000001 = 4L, subject_000002 = 4L))
  expect_equal(model$subject_000001@data@observed$Time, 0:3)
  expect_equal(model$subject_000002@data@observed$Time, 0:3)
  expect_false(anyNA(model$subject_000001@data@observed))
  expect_false(anyNA(model$subject_000002@data@observed))
})

test_that("multisubject likelihood is a subject-order-invariant sum", {
  build_fixed <- function(input) {
    model <- ctsgimme:::build_multisubject_ou_model(
      input,
      id_col = "id",
      drift_matrix = matrix("A_1,1", 1, 1),
      nvar = 1,
      varnames = "x",
      ME.var = matrix(0.1, 1, 1),
      PE.var = matrix(1, 1, 1),
      P0.values = matrix(1, 1, 1)
    )
    model$A$free[,] <- FALSE
    model$A$values[,] <- -0.4
    OpenMx::mxRun(model, silent = TRUE)
  }

  forward <- build_fixed(multisubject_test_data(FALSE))
  reverse <- build_fixed(multisubject_test_data(TRUE))
  child_objectives <- paste0(names(forward@submodels), ".fitfunction")
  child_sum <- sum(unlist(forward$output$algebras[child_objectives]))

  expect_equal(forward$output$fit, child_sum, tolerance = 1e-10)
  expect_equal(forward$output$fit, reverse$output$fit, tolerance = 1e-10)
})

test_that("explicit subgroup P0 values are validated", {
  input <- multisubject_test_data()
  build <- function(P0.values) {
    ctsgimme:::build_multisubject_ou_model(
      input,
      id_col = "id",
      drift_matrix = matrix("A_1,1", 1, 1),
      nvar = 1,
      varnames = "x",
      ME.var = matrix(0.1, 1, 1),
      PE.var = matrix(1, 1, 1),
      P0.values = P0.values
    )
  }

  expect_error(build(matrix(-1, 1, 1)), "positive definite", fixed = TRUE)
  expect_error(build(diag(2)), "one row and column", fixed = TRUE)
})

test_that("multisubject blocks reject unordered within-subject times", {
  input <- multisubject_test_data()
  input$Time[input$id == "B"] <- c(10, 12, 11, 13)

  expect_error(
    ctsgimme:::build_multisubject_ou_model(
      input,
      id_col = "id",
      drift_matrix = matrix("A_1,1", 1, 1),
      nvar = 1,
      varnames = "x",
      ME.var = matrix(0.1, 1, 1),
      PE.var = matrix(1, 1, 1),
      P0.values = matrix(1, 1, 1)
    ),
    "ordered by time",
    fixed = TRUE
  )
})

test_that("optional subgroup model construction failures are contained", {
  context <- list(
    id = "id",
    nvar = 1L,
    varnames = "x",
    ME.var = matrix(0.1, 1, 1),
    PE.var = matrix(1, 1, 1),
    ME.free = matrix(FALSE, 1, 1),
    PE.free = matrix(FALSE, 1, 1),
    time_col = "Time"
  )
  testthat::local_mocked_bindings(
    build_multisubject_ou_model = function(...) stop("intentional build failure"),
    .package = "ctsgimme"
  )

  expect_message(
    result <- ctsgimme:::.ctsgimme_fit_subgroup_model(
      context = context,
      new.data = multisubject_test_data(),
      DRIFT = matrix("A_1,1", 1, 1),
      G.DRIFT = matrix("0", 1, 1),
      subgroup = 1L,
      subgroup_dir = tempdir(),
      time.intervals = 1
    ),
    "intentional build failure",
    fixed = TRUE
  )
  expect_null(result)
})

test_that("subgroup transition intervals are validated", {
  expect_invisible(ctsgimme:::.ctsgimme_validate_time_intervals(c(0, 0.5, 2)))
  expect_error(
    ctsgimme:::.ctsgimme_validate_time_intervals(numeric()),
    "nonempty numeric vector",
    fixed = TRUE
  )
  expect_error(
    ctsgimme:::.ctsgimme_validate_time_intervals(c(1, -0.5)),
    "nonnegative",
    fixed = TRUE
  )
  expect_error(
    ctsgimme:::.ctsgimme_validate_time_intervals(c(1, Inf)),
    "finite",
    fixed = TRUE
  )
})

test_that("directed subgroup plot labels and colors follow plotted edges", {
  parameters <- data.frame(
    row = c(1L, 2L),
    col = c(2L, 1L),
    Estimate = c(0.4, -0.3),
    Std.Error = c(0.05, 0.05)
  )
  drift <- matrix("0", 2, 2)
  drift[1, 2] <- "A_1,2"
  drift[2, 1] <- "A_2,1"
  group_drift <- matrix("0", 2, 2)
  group_drift[2, 1] <- "A_2,1"

  components <- ctsgimme:::.ctsgimme_subgroup_plot_components(
    parameters,
    drift,
    group_drift,
    2L
  )
  colors <- matrix(components$edge.colors, 2, 2)

  expect_equal(components$effects[2, 1], 0.4)
  expect_equal(components$effects[1, 2], -0.3)
  expect_identical(components$edge.labels[2, 1], "* (0.4)")
  expect_identical(components$edge.labels[1, 2], "* (-0.3)")
  expect_identical(colors[2, 1], "blue")
  expect_identical(colors[1, 2], "gray")
})

test_that("real multisubject subgroup fit writes every requested artifact", {
  skip_if_not_installed("expm")
  make_series <- function(index, n = 45L) {
    innovations <- sin(seq_len(n) * (0.57 + index * 0.04)) +
      0.35 * cos(seq_len(n) * (1.11 + index * 0.03))
    x <- numeric(n)
    x[[1L]] <- c(0.25, -0.35, 0.15)[[index]]
    for (i in 2:n) {
      x[[i]] <- 0.58 * x[[i - 1L]] + 0.42 * innovations[[i]]
    }
    x
  }
  input <- do.call(rbind, lapply(seq_len(3L), function(index) {
    n <- 45L
    data.frame(
      id = LETTERS[[index]],
      x = make_series(index, n),
      Time = c(0, 100, 250)[[index]] + seq.int(0, n - 1L)
    )
  }))
  context <- list(
    id = "id",
    nvar = 1L,
    varnames = "x",
    ME.var = matrix(0.05, 1, 1),
    PE.var = matrix(0.5, 1, 1),
    ME.free = matrix(FALSE, 1, 1),
    PE.free = matrix(FALSE, 1, 1),
    time_col = "Time"
  )
  output_directory <- tempfile("ctsgimme-real-subgroup-")
  dir.create(output_directory)
  on.exit(unlink(output_directory, recursive = TRUE), add = TRUE)

  fit <- ctsgimme:::.ctsgimme_fit_subgroup_model(
    context = context,
    new.data = input,
    DRIFT = matrix("A_1,1", 1, 1),
    G.DRIFT = matrix("0", 1, 1),
    subgroup = 7L,
    subgroup_dir = output_directory,
    time.intervals = c(0.5, 2)
  )
  paths <- file.path(output_directory, c(
    "Subgroup 7 Params.png",
    "Subgroup 7 Delta_t = 0.5.png",
    "Subgroup 7 Delta_t = 2.png",
    "Subgroup_7Model.RDS"
  ))

  expect_true(all(file.exists(paths)))
  expect_true(all(file.info(paths)$size > 1000))
  png_signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
  for (path in paths[seq_len(3L)]) {
    expect_identical(readBin(path, "raw", 8), png_signature)
  }

  saved <- readRDS(paths[[4L]])
  expect_s4_class(fit, "MxModel")
  expect_s4_class(saved, "MxModel")
  expect_equal(saved$output$status$code, 0)
  expect_identical(attr(saved, "ctsgimme.subject.ids"), c("A", "B", "C"))
  expect_identical(
    attr(saved, "ctsgimme.subgroup.likelihood"),
    "multisubject"
  )
  expect_length(saved@submodels, 3L)
  expect_s4_class(saved@fitfunction, "MxFitFunctionMultigroup")
  groups <- sprintf("subject_%06d.fitfunction", seq_len(3L))
  expect_identical(saved@fitfunction@groups, groups)
  expect_true(all(c("A", "Q", "R", "x0", "P0") %in% names(saved@matrices)))
  expect_true(all(vapply(
    saved@submodels,
    function(model) model@data@observed$Time[[1L]] == 0,
    logical(1)
  )))
  expect_equal(
    saved$output$fit,
    sum(unlist(saved$output$algebras[groups])),
    tolerance = 1e-8
  )
  parameters <- subset(summary(saved)$parameters, matrix == "A")
  expect_equal(nrow(parameters), 1L)
  expect_true(is.finite(parameters$Estimate))
})
