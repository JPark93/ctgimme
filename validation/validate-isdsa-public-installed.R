args <- commandArgs(trailingOnly = TRUE)
named_args <- args[grepl("^--[^=]+=", args)]
options <- setNames(
  sub("^--[^=]+=", "", named_args),
  sub("^--([^=]+)=.*$", "\\1", named_args)
)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE))
} else {
  normalizePath("validation", mustWork = TRUE)
}
required <- "output"
missing <- setdiff(required, names(options))
if (length(missing)) {
  stop("Missing command argument(s): ", paste0("--", missing, collapse = ", "))
}

fixture_file <- normalizePath(
  if ("fixture" %in% names(options)) options[["fixture"]] else
    file.path(script_dir, "fixtures", "isdsa_jpark26_fixture.rds"),
  mustWork = TRUE
)
output_dir <- normalizePath(options[["output"]], mustWork = FALSE)
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
  stop("The output directory must be new or empty: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(ctgimme))
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required.")
if (!identical(as.character(packageVersion("ctgimme")), "0.1.0")) {
  stop("This acceptance test requires installed ctgimme 0.1.0.")
}
fixture_sha256 <- toupper(digest::digest(
  fixture_file,
  algo = "sha256",
  file = TRUE,
  serialize = FALSE
))
if (!identical(
  fixture_sha256,
  "5ADD83906470012F9DCE2B0F89098683B02370A42E698EDD018C9A8B24412985"
)) {
  stop("The ISDSA fixture SHA-256 does not match the frozen fixture.")
}
installed_package_path <- normalizePath(find.package("ctgimme"), mustWork = TRUE)
fixture <- readRDS(fixture_file)
if (nrow(fixture$sim) != 4000L || length(unique(fixture$sim$ID)) != 20L) {
  stop("The supplied fixture is not the expected 20-person ISDSA fixture.")
}

call_args <- fixture$ctsgimme_arguments
call_args$dataframe <- fixture$sim
call_args$cores <- 2L
call_args$directory <- output_dir
call_args$subgroup.method <- "pam"
call_args$max.subgroups <- 2L
call_args$subgroup.model <- TRUE
call_args$time.intervals <- c(0.25, 1, 2)
call_args$keep.intermediate <- FALSE
call_args$verbose <- FALSE

previous_threads <- OpenMx::mxOption(NULL, "Number of Threads")
OpenMx::mxOption(NULL, "Number of Threads", 1L)
on.exit(
  OpenMx::mxOption(NULL, "Number of Threads", previous_threads),
  add = TRUE
)
set.seed(20260819)
messages <- character()
warnings <- character()
started <- proc.time()[[3L]]
stdout <- utils::capture.output(
  result <- withCallingHandlers(
    do.call(ctgimme, call_args),
    message = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    },
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  ),
  type = "output"
)
elapsed <- proc.time()[[3L]] - started
if (length(stdout) || length(messages) || length(warnings)) {
  stop(
    "verbose = FALSE was not quiet. stdout=", length(stdout),
    ", messages=", length(messages), ", warnings=", length(warnings),
    if (length(warnings)) paste0("; ", paste(warnings, collapse = " | ")) else ""
  )
}

membership <- result$membership
subject_ids <- as.character(seq_along(fixture$truth$membership))
if (is.null(membership) || is.null(names(membership)) ||
    !all(subject_ids %in% names(membership))) {
  stop("The public result does not expose complete named membership.")
}
membership <- membership[subject_ids]
truth_membership <- stats::setNames(fixture$truth$membership, subject_ids)
if (!identical(
  outer(unname(membership), unname(membership), `==`),
  outer(unname(truth_membership), unname(truth_membership), `==`)
)) {
  stop("The public installed run did not recover the true ISDSA partition.")
}

contingency <- table(predicted = membership, truth = truth_membership)
mapping <- apply(contingency, 1L, function(counts) {
  as.integer(colnames(contingency)[[which.max(counts)]])
})
if (anyDuplicated(mapping) || any(apply(contingency, 1L, max) != 10L)) {
  stop("The recovered subgroup labels do not map one-to-one onto truth.")
}

group_structure <- readRDS(file.path(output_dir, "GStruc.RDS"))
if (!isTRUE(all.equal(
  group_structure != "0",
  fixture$truth$group_model != 0,
  check.attributes = FALSE
))) {
  stop("The public installed run did not recover the true group support.")
}

subgroup_rows <- list()
for (predicted_group in sort(unique(unname(membership)))) {
  true_group <- unname(mapping[[as.character(predicted_group)]])
  subgroup_structure <- readRDS(
    file.path(output_dir, paste0("SGStruc", predicted_group, ".RDS"))
  )
  if (!isTRUE(all.equal(
    subgroup_structure != "0",
    fixture$truth$subgroup_models[[true_group]] != 0,
    check.attributes = FALSE
  ))) {
    stop("A recovered subgroup structure differs from its true support.")
  }

  subgroup_dir <- file.path(
    output_dir,
    "Models",
    paste0("Subgroup ", predicted_group)
  )
  model_file <- file.path(
    subgroup_dir,
    paste0("Subgroup_", predicted_group, "Model.RDS")
  )
  subgroup_fit <- readRDS(model_file)
  expected_ids <- names(membership)[membership == predicted_group]
  if (subgroup_fit@output$status$code != 0 ||
      length(subgroup_fit@submodels) != 10L ||
      !methods::is(subgroup_fit@fitfunction, "MxFitFunctionMultigroup") ||
      length(OpenMx::omxGetParameters(subgroup_fit)) != 10L ||
      !identical(attr(subgroup_fit, "ctgimme.subgroup.likelihood"), "multisubject") ||
      !identical(attr(subgroup_fit, "ctgimme.subject.ids"), expected_ids)) {
    stop("A public subgroup RDS does not have the expected joint-model architecture.")
  }
  child_rows <- vapply(
    subgroup_fit@submodels,
    function(model) nrow(model@data@observed),
    integer(1)
  )
  child_starts <- vapply(
    subgroup_fit@submodels,
    function(model) model@data@observed$Time[[1L]],
    numeric(1)
  )
  groups <- subgroup_fit@fitfunction@groups
  if (any(child_rows != 200L) || any(child_starts != 0) ||
      abs(subgroup_fit@output$fit -
          sum(unlist(subgroup_fit@output$algebras[groups]))) > 1e-6) {
    stop("A public subgroup RDS has invalid child likelihood blocks.")
  }

  expected_plots <- file.path(subgroup_dir, c(
    paste0("Subgroup ", predicted_group, " Params.png"),
    paste0(
      "Subgroup ", predicted_group, " Delta_t = ",
      call_args$time.intervals,
      ".png"
    )
  ))
  if (!all(file.exists(c(model_file, expected_plots))) ||
      !all(file.info(c(model_file, expected_plots))$size > 1000)) {
    stop("A public subgroup output set is missing or empty.")
  }
  subgroup_rows[[length(subgroup_rows) + 1L]] <- data.frame(
    predicted_group = predicted_group,
    true_group = true_group,
    status_code = subgroup_fit@output$status$code,
    n_shared_free_parameters = length(OpenMx::omxGetParameters(subgroup_fit)),
    n_subject_likelihoods = length(subgroup_fit@submodels),
    stringsAsFactors = FALSE
  )
}

individual_dir <- file.path(output_dir, "Models", "Individuals")
individual_files <- file.path(
  individual_dir,
  paste0("FinalModel_", subject_ids, ".RDS")
)
individual_plots <- file.path(
  individual_dir,
  paste0("FinalModel_", subject_ids, ".PNG")
)
if (!all(file.exists(c(individual_files, individual_plots)))) {
  stop("The public installed run did not create 20 final individual outputs.")
}
individual_rows <- vector("list", length(individual_files))
for (index in seq_along(individual_files)) {
  fit <- readRDS(individual_files[[index]])
  estimated_support <- as.matrix(fit@matrices[["A"]]@free)
  true_support <- fixture$truth$individual_models[[index]] != 0
  historical_support <- true_support
  if (index == 5L) {
    # The retained installed 0.0.8 and 0.0.10 practice runs both selected
    # this finite-sample addition (X5 -> X4) for subject 5.
    historical_support[4L, 5L] <- TRUE
  }
  missing_true <- sum(!estimated_support & true_support)
  extra_false <- sum(estimated_support & !true_support)
  if (fit@output$status$code != 0) {
    stop("An individual model failed convergence: ", index)
  }
  individual_rows[[index]] <- data.frame(
    id = index,
    status_code = fit@output$status$code,
    n_estimated_paths = sum(estimated_support),
    n_true_paths = sum(true_support),
    missing_true_paths = missing_true,
    extra_false_paths = extra_false,
    exact_support = missing_true == 0L && extra_false == 0L,
    matches_installed_0_0_10_support = identical(
      unname(estimated_support),
      unname(historical_support)
    ),
    stringsAsFactors = FALSE
  )
}
individual_summary <- do.call(rbind, individual_rows)
if (!all(individual_summary$matches_installed_0_0_10_support)) {
  stop("Individual supports differ from the installed 0.0.10 practice baseline.")
}
if (sum(individual_summary$missing_true_paths) != 0L ||
    sum(individual_summary$extra_false_paths) > 1L ||
    sum(individual_summary$exact_support) < 19L) {
  stop("Individual support recovery fell outside the cross-version release bound.")
}

intermediate_models <- c(
  list.files(
    file.path(output_dir, "Models"),
    pattern = "^Model_.*\\.RDS$",
    full.names = TRUE
  ),
  unlist(lapply(
    list.dirs(file.path(output_dir, "Models"), recursive = FALSE),
    function(path) list.files(
      path,
      pattern = "^Model_.*\\.RDS$",
      full.names = TRUE
    )
  )),
  list.files(
    file.path(output_dir, "MIs"),
    pattern = "^MI_.*\\.RDS$",
    full.names = TRUE
  )
)
if (length(intermediate_models)) {
  stop("keep.intermediate = FALSE left shared-search model/MI artifacts.")
}

required_top_level <- file.path(output_dir, c(
  "GStruc.RDS",
  "SGStruc1.RDS",
  "SGStruc2.RDS",
  "subgroup_membership.csv",
  "subgroup_detection.rds",
  "Subgroups Plot.png"
))
if (!all(file.exists(required_top_level))) {
  stop("The public installed run is missing required top-level artifacts.")
}

subgroup_summary <- do.call(rbind, subgroup_rows)
utils::write.csv(
  subgroup_summary,
  file.path(output_dir, "public_acceptance_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  individual_summary,
  file.path(output_dir, "individual_support_summary.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    package_version = as.character(packageVersion("ctgimme")),
    package_path = installed_package_path,
    fixture = fixture_file,
    fixture_sha256 = fixture_sha256,
    elapsed_seconds = elapsed,
    membership = membership,
    truth_membership = truth_membership,
    label_mapping = mapping,
    subgroup_summary = subgroup_summary,
    individual_summary = individual_summary,
    session_info = utils::sessionInfo()
  ),
  file.path(output_dir, "public_acceptance_results.rds")
)
capture.output(utils::sessionInfo(), file = file.path(output_dir, "public_session-info.txt"))
writeLines(
  c(
    "# Installed ctgimme 0.1.0 public ISDSA acceptance",
    "",
    paste0("Elapsed seconds: ", format(elapsed, digits = 6)),
    "Membership, group support, and both subgroup supports matched fixture truth.",
    "All 20 individual support masks matched the installed 0.0.10 practice run.",
    paste0(
      sum(individual_summary$exact_support),
      "/20 individual supports were exact; all true paths were retained and ",
      sum(individual_summary$extra_false_paths),
      " extra path(s) appeared across all individuals."
    ),
    "Both subgroup RDS files contained one",
    "ten-subject multigroup likelihood and all requested PNG/RDS artifacts.",
    "The complete public run was silent with verbose = FALSE."
  ),
  file.path(output_dir, "PUBLIC_ACCEPTANCE.md"),
  useBytes = TRUE
)
message("Installed public ISDSA acceptance passed: ", output_dir)
