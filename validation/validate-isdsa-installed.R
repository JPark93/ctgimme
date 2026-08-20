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
reference_dir <- normalizePath(
  if ("reference" %in% names(options)) options[["reference"]] else
    file.path(script_dir, "reference-0.0.11"),
  mustWork = TRUE
)
output_dir <- normalizePath(
  options[["output"]],
  mustWork = FALSE
)
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
  stop("The output directory must be new or empty: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(ctgimme))
if (!requireNamespace("OpenMx", quietly = TRUE)) stop("OpenMx is required.")
if (!requireNamespace("expm", quietly = TRUE)) stop("expm is required.")
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required.")
if (!identical(as.character(packageVersion("ctgimme")), "0.0.12")) {
  stop("This acceptance test requires installed ctgimme 0.0.12.")
}
if (!identical(sort(getNamespaceExports("ctgimme")), c("ctgimme", "ctgimme_demo"))) {
  stop("The installed ctgimme exports are not the 0.0.12 public API.")
}

sha256_file <- function(path) {
  toupper(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}
sha256_text_lf <- function(path) {
  size <- file.info(path)$size
  bytes <- readBin(path, "raw", n = size)
  next_is_lf <- c(bytes[-1L] == as.raw(10L), FALSE)
  normalized <- bytes[!(bytes == as.raw(13L) & next_is_lf)]
  toupper(digest::digest(normalized, algo = "sha256", serialize = FALSE))
}
expected_hashes <- c(
  fixture = "5ADD83906470012F9DCE2B0F89098683B02370A42E698EDD018C9A8B24412985",
  fit_summary = "36A24EB4EDAE4AB8460C48BA2678A8EAD0EA6C597369C000B8FA402F3F51E3EF",
  parameter_estimates = "F1A7EF1BC628C7E1A3C874A96880E2BC76F26359FD2DCC10BCC26504064DB31F",
  subject_order = "29602E001B766270E476475434F66B23C59636C77521300BB821DF472609DDC1"
)
observed_hashes <- c(
  fixture = sha256_file(fixture_file),
  fit_summary = sha256_text_lf(file.path(reference_dir, "fit_summary.csv")),
  parameter_estimates = sha256_text_lf(file.path(reference_dir, "parameter_estimates.csv")),
  subject_order = sha256_text_lf(file.path(reference_dir, "subject_order_invariance.csv"))
)
if (!identical(observed_hashes, expected_hashes)) {
  stop("The ISDSA fixture or frozen 0.0.11 reference hash does not match.")
}
installed_package_path <- normalizePath(find.package("ctgimme"), mustWork = TRUE)

quick_membership <- unname(ctgimme_demo()$membership)
if (!identical(quick_membership, c(rep(1L, 4L), rep(2L, 4L)))) {
  stop("The installed quick demo did not recover its expected partition.")
}

reference_parameters <- utils::read.csv(
  file.path(reference_dir, "parameter_estimates.csv"),
  stringsAsFactors = FALSE
)
reference_summary <- utils::read.csv(
  file.path(reference_dir, "fit_summary.csv"),
  stringsAsFactors = FALSE
)
fixture <- readRDS(fixture_file)
if (nrow(fixture$sim) != 4000L || length(unique(fixture$sim$ID)) != 20L) {
  stop("The supplied fixture is not the expected 20-person ISDSA fixture.")
}
if (!identical(as.integer(table(fixture$truth$membership)), c(10L, 10L))) {
  stop("The supplied fixture does not contain two clone classes of ten.")
}
for (index in seq_along(fixture$truth$individual_models)) {
  subgroup <- fixture$truth$membership[[index]]
  if (!isTRUE(all.equal(
    fixture$truth$individual_models[[index]],
    fixture$truth$subgroup_models[[subgroup]],
    tolerance = 0
  ))) {
    stop("The fixture contains non-clone individual parameters.")
  }
}

previous_threads <- OpenMx::mxOption(NULL, "Number of Threads")
OpenMx::mxOption(NULL, "Number of Threads", 1L)
on.exit(
  OpenMx::mxOption(NULL, "Number of Threads", previous_threads),
  add = TRUE
)

varnames <- fixture$ctsgimme_arguments$varnames
nvar <- length(varnames)
ME.var <- fixture$ctsgimme_arguments$ME.var
PE.var <- fixture$ctsgimme_arguments$PE.var
time_intervals <- c(0.25, 1, 2)

drift_mask <- function(truth) {
  mask <- matrix("0", nrow(truth), ncol(truth))
  mask[truth != 0] <- paste0(
    "A_",
    row(truth)[truth != 0],
    ",",
    col(truth)[truth != 0]
  )
  mask
}

common_support <- Reduce(
  `&`,
  lapply(fixture$truth$subgroup_models, function(model) model != 0)
)
group_drift <- matrix("0", nvar, nvar)
group_drift[common_support] <- paste0(
  "A_",
  row(group_drift)[common_support],
  ",",
  col(group_drift)[common_support]
)

context <- list(
  id = "ID",
  nvar = nvar,
  varnames = varnames,
  ME.var = ME.var,
  PE.var = PE.var,
  ME.free = matrix(FALSE, nvar, nvar),
  PE.free = matrix(FALSE, nvar, nvar),
  time_col = "Time",
  verbose = FALSE
)

reorder_subjects <- function(data, reverse = FALSE) {
  ids <- unique(data$ID)
  if (reverse) ids <- rev(ids)
  result <- do.call(rbind, lapply(ids, function(id) {
    data[data$ID == id, , drop = FALSE]
  }))
  rownames(result) <- NULL
  result
}

extract_parameters <- function(fit, truth, subgroup, subject_order) {
  parameters <- subset(summary(fit)$parameters, matrix == "A")
  data.frame(
    subgroup = subgroup,
    subject_order = subject_order,
    row = as.integer(parameters$row),
    col = as.integer(parameters$col),
    parameter = paste0("A[", parameters$row, ",", parameters$col, "]"),
    truth = truth[cbind(as.integer(parameters$row), as.integer(parameters$col))],
    estimate = parameters$Estimate,
    std_error = parameters$Std.Error,
    stringsAsFactors = FALSE
  )
}

validate_architecture <- function(fit, expected_ids) {
  if (is.null(fit) || fit@output$status$code != 0) {
    stop("An installed-package ISDSA model did not converge with status 0.")
  }
  if (length(OpenMx::omxGetParameters(fit)) != 10L) {
    stop("An ISDSA model did not contain ten shared free parameters.")
  }
  if (length(fit@submodels) != 10L ||
      !methods::is(fit@fitfunction, "MxFitFunctionMultigroup")) {
    stop("An ISDSA model did not contain ten multigroup likelihood blocks.")
  }
  if (!identical(attr(fit, "ctgimme.subgroup.likelihood"), "multisubject") ||
      !identical(attr(fit, "ctgimme.subject.ids"), as.character(expected_ids))) {
    stop("Installed-model multisubject provenance is incorrect.")
  }
  child_rows <- vapply(
    fit@submodels,
    function(model) nrow(model@data@observed),
    integer(1)
  )
  child_starts <- vapply(
    fit@submodels,
    function(model) model@data@observed$Time[[1L]],
    numeric(1)
  )
  if (any(child_rows != 200L) || any(child_starts != 0)) {
    stop("Subject likelihood rows or local time origins are incorrect.")
  }
  groups <- fit@fitfunction@groups
  if (abs(fit@output$fit - sum(unlist(fit@output$algebras[groups]))) > 1e-6) {
    stop("The top-level likelihood is not the sum of subject likelihoods.")
  }
  invisible(TRUE)
}

results <- list()
fit_summaries <- list()
result_index <- 0L
started <- proc.time()[[3L]]
for (subgroup in seq_along(fixture$truth$subgroup_models)) {
  truth <- fixture$truth$subgroup_models[[subgroup]]
  mask <- drift_mask(truth)
  ids <- which(fixture$truth$membership == subgroup)
  subgroup_data <- fixture$sim[fixture$sim$ID %in% ids, , drop = FALSE]
  common_p0 <- ctgimme:::.ctgimme_initial_covariance(
    subgroup_data,
    varnames,
    ME.var
  )

  artifact_dir <- file.path(output_dir, paste0("Subgroup ", subgroup))
  dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
  fit_forward <- ctgimme:::.ctgimme_fit_subgroup_model(
    context = context,
    new.data = subgroup_data,
    DRIFT = mask,
    G.DRIFT = group_drift,
    subgroup = subgroup,
    subgroup_dir = artifact_dir,
    time.intervals = time_intervals
  )
  validate_architecture(fit_forward, ids)

  reverse_data <- reorder_subjects(subgroup_data, reverse = TRUE)
  reverse_model <- ctgimme:::build_multisubject_ou_model(
    data = reverse_data,
    id_col = "ID",
    drift_matrix = mask,
    nvar = nvar,
    varnames = varnames,
    ME.var = ME.var,
    PE.var = PE.var,
    ME.free = FALSE,
    PE.free = FALSE,
    model_name = "OUMod",
    lb = -10,
    ub = 10,
    time_col = "Time",
    P0.values = common_p0
  )
  fit_reverse <- ctgimme:::.ctgimme_mx_try_hard(
    reverse_model,
    verbose = FALSE,
    extraTries = 3
  )
  validate_architecture(fit_reverse, rev(ids))

  for (subject_order in c("forward", "reverse")) {
    fit <- if (subject_order == "forward") fit_forward else fit_reverse
    result_index <- result_index + 1L
    parameters <- extract_parameters(fit, truth, subgroup, subject_order)
    error <- parameters$estimate - parameters$truth
    results[[result_index]] <- parameters
    fit_summaries[[result_index]] <- data.frame(
      subgroup = subgroup,
      subject_order = subject_order,
      status_code = fit@output$status$code,
      minus2_log_likelihood = fit@output$fit,
      n_free_parameters = length(OpenMx::omxGetParameters(fit)),
      n_subject_likelihoods = length(fit@submodels),
      bias = mean(error),
      mae = mean(abs(error)),
      rmse = sqrt(mean(error^2)),
      max_abs_error = max(abs(error)),
      stringsAsFactors = FALSE
    )
  }

  expected_artifacts <- file.path(artifact_dir, c(
    paste0("Subgroup ", subgroup, " Params.png"),
    paste0("Subgroup ", subgroup, " Delta_t = ", time_intervals, ".png"),
    paste0("Subgroup_", subgroup, "Model.RDS")
  ))
  if (!all(file.exists(expected_artifacts)) ||
      !all(file.info(expected_artifacts)$size > 1000)) {
    stop("The installed subgroup fitter did not create all nonempty artifacts.")
  }
  png_signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
  for (path in expected_artifacts[grepl("\\.png$", expected_artifacts)]) {
    if (!identical(readBin(path, "raw", 8L), png_signature)) {
      stop("An installed subgroup plot is not a valid PNG: ", path)
    }
  }
}

parameters <- do.call(rbind, results)
fit_summary <- do.call(rbind, fit_summaries)
parameters$error <- parameters$estimate - parameters$truth
parameters$abs_error <- abs(parameters$error)

keys <- c("subgroup", "subject_order", "row", "col", "parameter")
comparison <- merge(
  parameters,
  reference_parameters[, c(keys, "estimate", "std_error")],
  by = keys,
  suffixes = c("_current", "_reference")
)
if (nrow(comparison) != nrow(parameters)) {
  stop("Current and reference parameter sets do not have identical keys.")
}
max_parameter_difference <- max(abs(
  comparison$estimate_current - comparison$estimate_reference
))
if (max_parameter_difference > 1e-5) {
  stop("Parameter parity exceeded 1e-5: ", max_parameter_difference)
}
max_standard_error_difference <- max(abs(
  comparison$std_error_current - comparison$std_error_reference
))
if (max_standard_error_difference > 1e-4) {
  stop("Standard-error parity exceeded 1e-4: ", max_standard_error_difference)
}

summary_comparison <- merge(
  fit_summary,
  reference_summary[, c("subgroup", "subject_order", "minus2_log_likelihood")],
  by = c("subgroup", "subject_order"),
  suffixes = c("_current", "_reference")
)
max_likelihood_difference <- max(abs(
  summary_comparison$minus2_log_likelihood_current -
    summary_comparison$minus2_log_likelihood_reference
))
if (max_likelihood_difference > 1e-4) {
  stop("Likelihood parity exceeded 1e-4: ", max_likelihood_difference)
}

forward <- parameters[parameters$subject_order == "forward", ]
reverse <- parameters[parameters$subject_order == "reverse", ]
order_comparison <- merge(
  forward[, c("subgroup", "row", "col", "estimate")],
  reverse[, c("subgroup", "row", "col", "estimate")],
  by = c("subgroup", "row", "col"),
  suffixes = c("_forward", "_reverse")
)
order_summary <- aggregate(
  abs(order_comparison$estimate_forward - order_comparison$estimate_reverse),
  by = order_comparison["subgroup"],
  FUN = max
)
names(order_summary)[[2L]] <- "max_abs_order_difference"
if (any(order_summary$max_abs_order_difference > 1e-6)) {
  stop("Subject-order estimate invariance exceeded 1e-6.")
}

forward_summary <- fit_summary[fit_summary$subject_order == "forward", ]
sg1 <- forward_summary[forward_summary$subgroup == 1L, ]
sg2 <- forward_summary[forward_summary$subgroup == 2L, ]
combined_error <- forward$error
truth_ok <-
  sg1$mae <= 0.04 && sg1$rmse <= 0.05 && sg1$max_abs_error <= 0.10 &&
  sg2$mae <= 0.03 && sg2$rmse <= 0.035 && sg2$max_abs_error <= 0.07 &&
  abs(mean(combined_error)) <= 0.01 &&
  mean(abs(combined_error)) <= 0.03 &&
  sqrt(mean(combined_error^2)) <= 0.04
if (!truth_ok) stop("Installed-package truth-recovery bounds were not met.")

elapsed <- proc.time()[[3L]] - started
utils::write.csv(fit_summary, file.path(output_dir, "fit_summary.csv"), row.names = FALSE)
utils::write.csv(parameters, file.path(output_dir, "parameter_estimates.csv"), row.names = FALSE)
utils::write.csv(order_summary, file.path(output_dir, "subject_order_invariance.csv"), row.names = FALSE)
saveRDS(
  list(
    package_version = as.character(packageVersion("ctgimme")),
    package_path = installed_package_path,
    fixture = fixture_file,
    fixture_sha256 = observed_hashes[["fixture"]],
    reference = reference_dir,
    reference_sha256 = observed_hashes[c("fit_summary", "parameter_estimates", "subject_order")],
    elapsed_seconds = elapsed,
    max_parameter_difference = max_parameter_difference,
    max_standard_error_difference = max_standard_error_difference,
    max_likelihood_difference = max_likelihood_difference,
    fit_summary = fit_summary,
    parameter_estimates = parameters,
    subject_order_invariance = order_summary,
    session_info = utils::sessionInfo()
  ),
  file.path(output_dir, "validation_results.rds")
)
capture.output(utils::sessionInfo(), file = file.path(output_dir, "session-info.txt"))

report <- c(
  "# Installed ctgimme 0.0.12 ISDSA acceptance",
  "",
  paste0("Installed package: ", installed_package_path),
  paste0("Fixture SHA-256: ", observed_hashes[["fixture"]]),
  paste0("Elapsed seconds: ", format(elapsed, digits = 6)),
  paste0("Maximum estimate difference from 0.0.11 reference: ", format(max_parameter_difference, digits = 12)),
  paste0("Maximum standard-error difference from 0.0.11 reference: ", format(max_standard_error_difference, digits = 12)),
  paste0("Maximum -2LL difference from 0.0.11 reference: ", format(max_likelihood_difference, digits = 12)),
  paste0("Maximum forward/reverse estimate difference: ", format(max(order_summary$max_abs_order_difference), digits = 12)),
  "",
  "All convergence, architecture, parity, order-invariance, truth-recovery,",
  "quick-demo, and PNG/RDS artifact gates passed."
)
writeLines(report, file.path(output_dir, "REPORT.md"), useBytes = TRUE)
message("Installed-package ISDSA acceptance passed: ", output_dir)
