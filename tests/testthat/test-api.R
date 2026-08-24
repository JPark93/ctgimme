test_that("the public API exposes documented defaults", {
  api <- formals(ctgimme)

  expect_identical(
    tail(names(api), 4L),
    c("scale.data", "ME.free", "PE.free", "verbose")
  )
  expect_identical(eval(api$cores), 1)
  expect_identical(eval(api$scale.data), FALSE)
  expect_identical(eval(api$ME.free), FALSE)
  expect_identical(eval(api$PE.free), FALSE)
  expect_identical(eval(api$verbose), TRUE)
  expect_identical(eval(api$subgroup.method), c("pam", "legacy"))
  expect_false(any(c(
    "subgroup.model.method",
    "subgroup.time.mode",
    "measurement.schedule",
    "cycle.interval",
    "insert.na.rows",
    "subject.gap"
  ) %in% names(api)))
})

test_that("the renamed namespace exposes only the ctgimme public API", {
  exports <- getNamespaceExports("ctgimme")

  expect_true("ctgimme" %in% exports)
  expect_true("ctgimme_demo" %in% exports)
  expect_false("ctsgimme" %in% exports)
  expect_false("ctsgimme_demo" %in% exports)
})

test_that("parallel requests are bounded only by the subject count", {
  expect_identical(ctgimme:::.ctgimme_resolve_cores(1, 20), 1L)
  expect_identical(ctgimme:::.ctgimme_resolve_cores(2, 20), 2L)
  expect_identical(ctgimme:::.ctgimme_resolve_cores(6, 20), 6L)
  expect_identical(ctgimme:::.ctgimme_resolve_cores(20, 20), 20L)
  expect_message(
    expect_identical(ctgimme:::.ctgimme_resolve_cores(6, 1), 1L),
    "number of subjects",
    fixed = TRUE
  )
  expect_no_message(
    resolved <- ctgimme:::.ctgimme_resolve_cores(6, 1, verbose = FALSE)
  )
  expect_identical(resolved, 1L)
  expect_error(
    ctgimme:::.ctgimme_resolve_cores(0, 20),
    "one positive integer",
    fixed = TRUE
  )
})

test_that("core-adjustment messages honor verbose", {
  expect_no_message(
    resolved <- ctgimme:::.ctgimme_resolve_cores(6, 4, verbose = FALSE)
  )
  expect_identical(resolved, 4L)
})

test_that("large explicit core requests have no package-wide maximum", {
  expect_identical(ctgimme:::.ctgimme_resolve_cores(64, 100), 64L)
})

test_that("worker-cluster creation preserves resolved counts above two", {
  tracker <- new.env(parent = emptyenv())
  tracker$cores <- NULL
  fake_cluster <- structure(
    list(id = "four-worker-pool"),
    class = "mock_cluster"
  )

  cluster <- testthat::with_mocked_bindings(
    testthat::with_mocked_bindings(
      ctgimme:::.ctgimme_make_worker_cluster(4L),
      clusterCall = function(...) list(TRUE),
      .package = "parallel"
    ),
    makeCluster = function(cores, type) {
      tracker$cores <- cores
      expect_identical(type, "PSOCK")
      fake_cluster
    },
    clusterEvalQ = function(cluster, expression) rep(TRUE, tracker$cores),
    .ctgimme_register_worker_cluster = function(cluster) cluster,
    .ctgimme_stop_worker_cluster = function(cluster) invisible(TRUE),
    .package = "ctgimme"
  )

  expect_identical(tracker$cores, 4L)
  expect_identical(cluster, fake_cluster)
})

test_that("the quick demo deterministically recovers two subgroups", {
  result <- ctgimme_demo()

  expect_identical(length(result$membership), 8L)
  expect_identical(unname(result$membership[1:4]), rep(1L, 4))
  expect_identical(unname(result$membership[5:8]), rep(2L, 4))
  expect_identical(result$selection$selected, c(TRUE, TRUE, FALSE))
  expect_identical(result$candidates$k, 2L)
})

test_that("basic invalid inputs fail before estimation", {
  expect_error(ctgimme(), "varnames must be supplied", fixed = TRUE)
  expect_error(
    ctgimme(verbose = NA),
    "verbose must be supplied as TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    ctgimme(verbose = 1),
    "verbose must be supplied as TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    ctgimme(verbose = c(TRUE, FALSE)),
    "verbose must be supplied as TRUE or FALSE",
    fixed = TRUE
  )

  input <- data.frame(id = 1, time = 0, x = 0)
  expect_error(
    ctgimme(
      varnames = "x",
      dataframe = input,
      id = "id",
      time = "time",
      directory = tempfile(),
      scale.data = "FALSE"
    ),
    "scale.data must be supplied as TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    ctgimme(
      varnames = "x",
      dataframe = input,
      id = "id",
      time = "time",
      directory = tempfile(),
      sub.sig.thrsh = 0.55,
      max.subgroups = NULL
    ),
    "max.subgroups must be supplied",
    fixed = TRUE
  )
})

test_that("documented data and control constraints are enforced", {
  valid <- data.frame(
    id = rep(c("S1", "S2"), each = 2L),
    time = rep(0:1, 2L),
    x = c(1, NA, 2, 3)
  )
  expect_invisible(
    ctgimme:::.ctgimme_validate_inputs(
      "x", valid, "id", "time", tempfile()
    )
  )
  expect_error(
    ctgimme:::.ctgimme_validate_inputs(
      c("x", "x"), valid, "id", "time", tempfile()
    ),
    "distinct column names",
    fixed = TRUE
  )
  expect_error(
    ctgimme:::.ctgimme_validate_inputs(
      "x", transform(valid, id = ".."), "id", "time", tempfile()
    ),
    "safe for use in filenames",
    fixed = TRUE
  )
  case_collision <- transform(
    valid,
    id = rep(c("SubjectA", "subjecta"), each = 2L)
  )
  expect_error(
    ctgimme:::.ctgimme_validate_inputs(
      "x", case_collision, "id", "time", tempfile()
    ),
    "distinct when compared without case",
    fixed = TRUE
  )
  unordered <- valid[c(2L, 1L, 3L, 4L), ]
  expect_error(
    ctgimme:::.ctgimme_validate_inputs(
      "x", unordered, "id", "time", tempfile()
    ),
    "ordered by time",
    fixed = TRUE
  )
  expect_error(
    ctgimme:::.ctgimme_validate_inputs(
      "x", transform(valid, x = as.character(x)), "id", "time", tempfile()
    ),
    "must be numeric",
    fixed = TRUE
  )
  expect_error(
    ctgimme:::.ctgimme_validate_controls(
      1.1, 1, 0.05, TRUE, 0.05, 0.01, FALSE, TRUE
    ),
    "sig.thrsh must be one finite numeric value in [0, 1]",
    fixed = TRUE
  )
  expect_error(
    ctgimme:::.ctgimme_validate_controls(
      0.55, 1, 0.05, TRUE, 0.05, 0.01, "FALSE", TRUE
    ),
    "keep.intermediate must be supplied as TRUE or FALSE",
    fixed = TRUE
  )
})

test_that("verbose suppresses package progress without changing results", {
  subject_ids <- c("S1", "S2")
  membership <- stats::setNames(c(1L, 1L), subject_ids)
  output_directory <- tempfile("ctgimme-verbose-")

  testthat::local_mocked_bindings(
    .ctgimme_validate_inputs = function(...) invisible(TRUE),
    .ctgimme_validate_subgroup_options = function(...) invisible(TRUE),
    .ctgimme_load_dependencies = function(...) invisible(TRUE),
    .ctgimme_prepare_context = function(...) {
      arguments <- list(...)
      list(
        directory = output_directory,
        ids = subject_ids,
        cores = 1L,
        verbose = arguments[[length(arguments)]]
      )
    },
    .ctgimme_run_group_stage = function(...) matrix("0", 1, 1),
    .ctgimme_detect_subgroups = function(...) {
      list(
        membership = membership,
        clustering = NULL,
        method = "pam",
        diagnostics = list()
      )
    },
    .ctgimme_write_subgroup_membership = function(...) list(),
    .ctgimme_save_rds = function(...) invisible(NULL),
    .ctgimme_run_subgroup_stages = function(...) invisible(NULL),
    .ctgimme_cleanup = function(...) invisible(NULL),
    .package = "ctgimme"
  )

  run <- function(verbose) {
    ctgimme(
      varnames = "x",
      dataframe = data.frame(id = subject_ids, time = 0, x = 0),
      id = "id",
      time = "time",
      directory = output_directory,
      verbose = verbose
    )
  }

  expect_silent(quiet <- run(FALSE))
  messages <- character()
  loud <- withCallingHandlers(
    run(TRUE),
    message = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl("Assigned 2 subject", messages, fixed = TRUE)))
  expect_true(any(grepl("Complete", messages, fixed = TRUE)))
  expect_identical(quiet$membership, loud$membership)
})

test_that("a quiet two-worker public run exports every fitter dependency", {
  make_series <- function(offset, n = 12L) {
    values <- numeric(n)
    values[[1L]] <- 0.2 * offset
    for (index in 2:n) {
      values[[index]] <-
        0.55 * values[[index - 1L]] + sin((index + offset) * 0.7)
    }
    values
  }
  n_time <- 12L
  input <- data.frame(
    id = rep(c("P1", "P2"), each = n_time),
    time = rep(0:(n_time - 1L), 2L),
    x = c(make_series(1), make_series(2))
  )
  output_directory <- tempfile("ctgimme-real-parallel-")
  on.exit(unlink(output_directory, recursive = TRUE), add = TRUE)

  expect_silent(
    result <- ctgimme(
      varnames = "x",
      dataframe = input,
      id = "id",
      time = "time",
      cores = 2,
      directory = output_directory,
      ME.var = 0.05,
      PE.var = 1,
      verbose = FALSE
    )
  )

  expect_identical(names(result$membership), c("P1", "P2"))
  expect_length(result$membership, 2L)
})

test_that("the public workflow records subgroup membership artifacts", {
  output_directory <- tempfile("ctgimme-public-membership-")
  subject_ids <- as.character(seq_len(6))
  membership <- stats::setNames(c(1L, 1L, 1L, 2L, 2L, 2L), subject_ids)
  subject_points <- rbind(
    c(0, 0), c(0, 0.1), c(0.1, 0),
    c(4, 4), c(4, 4.1), c(4.1, 4)
  )
  rownames(subject_points) <- subject_ids
  subject_distance <- as.matrix(stats::dist(subject_points, method = "manhattan"))
  pam_fit <- cluster::pam(subject_distance, k = 2, diss = TRUE)

  testthat::local_mocked_bindings(
    .ctgimme_validate_inputs = function(...) invisible(TRUE),
    .ctgimme_validate_subgroup_options = function(...) invisible(TRUE),
    .ctgimme_load_dependencies = function(...) invisible(TRUE),
    .ctgimme_prepare_context = function(...) {
      list(directory = output_directory, ids = subject_ids, cores = 1L)
    },
    .ctgimme_run_group_stage = function(...) matrix("0", 1, 1),
    .ctgimme_detect_subgroups = function(...) {
      list(
        membership = membership,
        clustering = pam_fit,
        method = "pam",
        diagnostics = list(distance = subject_distance)
      )
    },
    .ctgimme_run_subgroup_stages = function(...) invisible(NULL),
    .ctgimme_cleanup = function(...) invisible(NULL),
    .package = "ctgimme"
  )

  result <- ctgimme(
    varnames = "x",
    dataframe = data.frame(id = subject_ids, time = 0, x = 0),
    id = "id",
    time = "time",
    directory = output_directory,
    sub.sig.thrsh = 0.55,
    max.subgroups = 2
  )
  detection <- readRDS(file.path(output_directory, "subgroup_detection.rds"))
  artifacts <- attr(result, "ctgimme.membership.artifacts")

  expect_true(file.exists(file.path(output_directory, "Subgroups Plot.png")))
  expect_true(file.exists(file.path(output_directory, "subgroup_membership.csv")))
  expect_identical(attr(result, "ctgimme.membership"), membership)
  expect_identical(artifacts$plot.type, "distance-map")
  expect_identical(detection$membership, membership)
  expect_identical(detection$diagnostics$membership.artifacts, artifacts)
})

test_that("the public workflow reuses and stops one worker cluster", {
  tracker <- new.env(parent = emptyenv())
  tracker$created <- 0L
  tracker$stopped <- 0L
  tracker$group.cluster <- NULL
  tracker$subgroup.cluster <- NULL
  worker_cluster <- structure(list(id = "persistent"), class = "mock_cluster")
  subject_ids <- as.character(seq_len(4L))
  membership <- stats::setNames(rep(1L, 4L), subject_ids)

  testthat::local_mocked_bindings(
    .ctgimme_validate_inputs = function(...) invisible(TRUE),
    .ctgimme_validate_subgroup_options = function(...) invisible(TRUE),
    .ctgimme_load_dependencies = function(...) invisible(TRUE),
    .ctgimme_prepare_context = function(...) {
      list(
        directory = tempfile("ctgimme-pool-test-"),
        ids = subject_ids,
        cores = 4L
      )
    },
    .ctgimme_make_worker_cluster = function(cores) {
      tracker$created <- tracker$created + 1L
      expect_identical(cores, 4L)
      worker_cluster
    },
    stopCluster = function(cluster) {
      expect_identical(cluster, worker_cluster)
      tracker$stopped <- tracker$stopped + 1L
      invisible(NULL)
    },
    .ctgimme_run_group_stage = function(
        context, Galpha, sig.thrsh, worker_cluster = NULL) {
      tracker$group.cluster <- worker_cluster
      matrix("0", 1, 1)
    },
    .ctgimme_detect_subgroups = function(...) {
      list(
        membership = membership,
        clustering = NULL,
        method = "pam",
        diagnostics = list()
      )
    },
    .ctgimme_write_subgroup_membership = function(...) list(),
    .ctgimme_save_rds = function(...) invisible(NULL),
    .ctgimme_run_subgroup_stages = function(..., worker_cluster = NULL) {
      tracker$subgroup.cluster <- worker_cluster
      invisible(NULL)
    },
    .ctgimme_cleanup = function(...) invisible(NULL),
    .package = "ctgimme"
  )

  invisible(ctgimme(
    varnames = "x",
    dataframe = data.frame(id = subject_ids, time = 0, x = 0),
    id = "id",
    time = "time",
    directory = tempfile("ctgimme-pool-output-"),
    cores = 4,
    sub.sig.thrsh = 1
  ))

  expect_identical(tracker$created, 1L)
  expect_identical(tracker$stopped, 1L)
  expect_identical(tracker$group.cluster, worker_cluster)
  expect_identical(tracker$subgroup.cluster, worker_cluster)
})

test_that("the public workflow stops its worker cluster after an error", {
  tracker <- new.env(parent = emptyenv())
  tracker$created <- 0L
  tracker$stopped <- 0L
  worker_cluster <- structure(list(id = "error-cleanup"), class = "mock_cluster")

  testthat::local_mocked_bindings(
    .ctgimme_validate_inputs = function(...) invisible(TRUE),
    .ctgimme_validate_subgroup_options = function(...) invisible(TRUE),
    .ctgimme_load_dependencies = function(...) invisible(TRUE),
    .ctgimme_prepare_context = function(...) {
      list(
        directory = tempfile("ctgimme-error-pool-test-"),
        ids = c("1", "2"),
        cores = 2L
      )
    },
    .ctgimme_make_worker_cluster = function(cores) {
      tracker$created <- tracker$created + 1L
      worker_cluster
    },
    stopCluster = function(cluster) {
      expect_identical(cluster, worker_cluster)
      tracker$stopped <- tracker$stopped + 1L
      invisible(NULL)
    },
    .ctgimme_run_group_stage = function(...) {
      stop("intentional stage failure")
    },
    .package = "ctgimme"
  )

  expect_error(
    ctgimme(
      varnames = "x",
      dataframe = data.frame(id = c("1", "2"), time = 0, x = 0),
      id = "id",
      time = "time",
      directory = tempfile("ctgimme-error-pool-output-"),
      cores = 2,
      sub.sig.thrsh = 1
    ),
    "intentional stage failure"
  )

  expect_identical(tracker$created, 1L)
  expect_identical(tracker$stopped, 1L)
})
