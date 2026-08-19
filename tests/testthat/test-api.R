test_that("the public API exposes documented defaults", {
  api <- formals(ctsgimme)

  expect_identical(
    tail(names(api), 3L),
    c("scale.data", "ME.free", "PE.free")
  )
  expect_identical(eval(api$cores), 1)
  expect_identical(eval(api$scale.data), FALSE)
  expect_identical(eval(api$ME.free), FALSE)
  expect_identical(eval(api$PE.free), FALSE)
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

test_that("core requests are not unnecessarily capped during ordinary use", {
  check_setting <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA_character_)
  on.exit({
    if (is.na(check_setting)) {
      Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
    } else {
      do.call(Sys.setenv, setNames(list(check_setting), "_R_CHECK_LIMIT_CORES_"))
    }
  }, add = TRUE)
  Sys.unsetenv("_R_CHECK_LIMIT_CORES_")

  expect_identical(ctsgimme:::.ctsgimme_resolve_cores(1, 20), 1L)
  expect_identical(ctsgimme:::.ctsgimme_resolve_cores(2, 20), 2L)
  expect_identical(ctsgimme:::.ctsgimme_resolve_cores(6, 20), 6L)
  expect_message(
    expect_identical(ctsgimme:::.ctsgimme_resolve_cores(6, 1), 1L),
    "number of subjects",
    fixed = TRUE
  )
  expect_error(
    ctsgimme:::.ctsgimme_resolve_cores(0, 20),
    "one positive integer",
    fixed = TRUE
  )
})

test_that("CRAN-style checks limit execution to two workers", {
  check_setting <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA_character_)
  on.exit({
    if (is.na(check_setting)) {
      Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
    } else {
      do.call(Sys.setenv, setNames(list(check_setting), "_R_CHECK_LIMIT_CORES_"))
    }
  }, add = TRUE)

  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")
  expect_message(
    expect_identical(ctsgimme:::.ctsgimme_resolve_cores(6, 20), 2L),
    "CRAN-style package checking",
    fixed = TRUE
  )

  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "warn")
  expect_message(
    expect_identical(ctsgimme:::.ctsgimme_resolve_cores(6, 20), 2L),
    "CRAN-style package checking",
    fixed = TRUE
  )

  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "false")
  expect_identical(ctsgimme:::.ctsgimme_resolve_cores(6, 20), 6L)
})

test_that("the quick demo deterministically recovers two subgroups", {
  result <- ctsgimme_demo()

  expect_identical(length(result$membership), 8L)
  expect_identical(unname(result$membership[1:4]), rep(1L, 4))
  expect_identical(unname(result$membership[5:8]), rep(2L, 4))
  expect_identical(result$selection$selected, c(TRUE, TRUE, FALSE))
  expect_identical(result$candidates$k, 2L)
})

test_that("basic invalid inputs fail before estimation", {
  expect_error(ctsgimme(), "varnames must be supplied", fixed = TRUE)

  input <- data.frame(id = 1, time = 0, x = 0)
  expect_error(
    ctsgimme(
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
    ctsgimme(
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

test_that("the public workflow records subgroup membership artifacts", {
  output_directory <- tempfile("ctsgimme-public-membership-")
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
    .ctsgimme_validate_inputs = function(...) invisible(TRUE),
    .ctsgimme_validate_subgroup_options = function(...) invisible(TRUE),
    .ctsgimme_load_dependencies = function(...) invisible(TRUE),
    .ctsgimme_prepare_context = function(...) {
      list(directory = output_directory, ids = subject_ids, cores = 1L)
    },
    .ctsgimme_run_group_stage = function(...) matrix("0", 1, 1),
    .ctsgimme_detect_subgroups = function(...) {
      list(
        membership = membership,
        clustering = pam_fit,
        method = "pam",
        diagnostics = list(distance = subject_distance)
      )
    },
    .ctsgimme_run_subgroup_stages = function(...) invisible(NULL),
    .ctsgimme_cleanup = function(...) invisible(NULL),
    .package = "ctsgimme"
  )

  result <- ctsgimme(
    varnames = "x",
    dataframe = data.frame(id = subject_ids, time = 0, x = 0),
    id = "id",
    time = "time",
    directory = output_directory,
    sub.sig.thrsh = 0.55,
    max.subgroups = 2
  )
  detection <- readRDS(file.path(output_directory, "subgroup_detection.rds"))
  artifacts <- attr(result, "ctsgimme.membership.artifacts")

  expect_true(file.exists(file.path(output_directory, "Subgroups Plot.png")))
  expect_true(file.exists(file.path(output_directory, "subgroup_membership.csv")))
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
  subject_ids <- c("1", "2")
  membership <- stats::setNames(c(1L, 1L), subject_ids)

  testthat::local_mocked_bindings(
    .ctsgimme_validate_inputs = function(...) invisible(TRUE),
    .ctsgimme_validate_subgroup_options = function(...) invisible(TRUE),
    .ctsgimme_load_dependencies = function(...) invisible(TRUE),
    .ctsgimme_prepare_context = function(...) {
      list(directory = tempfile("ctsgimme-pool-test-"), ids = subject_ids, cores = 2L)
    },
    .ctsgimme_make_worker_cluster = function(cores) {
      tracker$created <- tracker$created + 1L
      expect_identical(cores, 2L)
      worker_cluster
    },
    stopCluster = function(cluster) {
      expect_identical(cluster, worker_cluster)
      tracker$stopped <- tracker$stopped + 1L
      invisible(NULL)
    },
    .ctsgimme_run_group_stage = function(
        context, Galpha, sig.thrsh, worker_cluster = NULL) {
      tracker$group.cluster <- worker_cluster
      matrix("0", 1, 1)
    },
    .ctsgimme_detect_subgroups = function(...) {
      list(
        membership = membership,
        clustering = NULL,
        method = "pam",
        diagnostics = list()
      )
    },
    .ctsgimme_write_subgroup_membership = function(...) list(),
    .ctsgimme_save_rds = function(...) invisible(NULL),
    .ctsgimme_run_subgroup_stages = function(..., worker_cluster = NULL) {
      tracker$subgroup.cluster <- worker_cluster
      invisible(NULL)
    },
    .ctsgimme_cleanup = function(...) invisible(NULL),
    .package = "ctsgimme"
  )

  invisible(ctsgimme(
    varnames = "x",
    dataframe = data.frame(id = subject_ids, time = 0, x = 0),
    id = "id",
    time = "time",
    directory = tempfile("ctsgimme-pool-output-"),
    cores = 2,
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
    .ctsgimme_validate_inputs = function(...) invisible(TRUE),
    .ctsgimme_validate_subgroup_options = function(...) invisible(TRUE),
    .ctsgimme_load_dependencies = function(...) invisible(TRUE),
    .ctsgimme_prepare_context = function(...) {
      list(
        directory = tempfile("ctsgimme-error-pool-test-"),
        ids = c("1", "2"),
        cores = 2L
      )
    },
    .ctsgimme_make_worker_cluster = function(cores) {
      tracker$created <- tracker$created + 1L
      worker_cluster
    },
    stopCluster = function(cluster) {
      expect_identical(cluster, worker_cluster)
      tracker$stopped <- tracker$stopped + 1L
      invisible(NULL)
    },
    .ctsgimme_run_group_stage = function(...) {
      stop("intentional stage failure")
    },
    .package = "ctsgimme"
  )

  expect_error(
    ctsgimme(
      varnames = "x",
      dataframe = data.frame(id = c("1", "2"), time = 0, x = 0),
      id = "id",
      time = "time",
      directory = tempfile("ctsgimme-error-pool-output-"),
      cores = 2,
      sub.sig.thrsh = 1
    ),
    "intentional stage failure"
  )

  expect_identical(tracker$created, 1L)
  expect_identical(tracker$stopped, 1L)
})
