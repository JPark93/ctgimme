.ctsgimme_extract_subgroup_score_features <- function(context, group_drift) {
  ids <- as.character(context$ids)
  features <- matrix(
    0,
    nrow = length(ids),
    ncol = length(context$param_names),
    dimnames = list(ids, context$param_names)
  )
  kept.ids <- character(0)

  for (subject_id in ids) {
    model_file <- file.path(
      context$directory,
      "Models",
      paste0("Model_", subject_id, ".RDS")
    )
    mi_file <- file.path(
      context$directory,
      "MIs",
      paste0("MI_", subject_id, ".RDS")
    )
    fit <- tryCatch(readRDS(model_file), error = function(e) NULL)
    mi <- tryCatch(readRDS(mi_file), error = function(e) NULL)
    if (is.null(fit) || is.null(mi) ||
        is.null(mi$MI.Full) || is.null(mi$EPC)) next

    estimates <- tryCatch(
      subset(summary(fit)$parameters, matrix == "A"),
      error = function(e) NULL
    )
    if (is.null(estimates)) next
    if (nrow(estimates)) {
      estimate_names <- if ("name" %in% names(estimates)) {
        as.character(estimates$name)
      } else {
        sprintf("OUMod.A[%s,%s]", estimates$row, estimates$col)
      }
      estimate_z <- estimates$Estimate / estimates$Std.Error
      estimate_z[!is.finite(estimate_z)] <- 0
      matched <- match(estimate_names, colnames(features))
      keep <- !is.na(matched)
      features[subject_id, matched[keep]] <- estimate_z[keep]
    }

    mi_values <- suppressWarnings(as.numeric(mi$MI.Full))
    epc_values <- suppressWarnings(as.numeric(mi$EPC))
    if (length(mi_values) != length(epc_values)) next
    score <- sign(epc_values) * sqrt(pmax(abs(mi_values), 0))
    score[!is.finite(score)] <- 0
    matched <- match(names(mi$MI.Full), colnames(features))
    keep <- !is.na(matched)
    features[subject_id, matched[keep]] <- score[keep]
    kept.ids <- c(kept.ids, subject_id)
  }

  cells <- t(vapply(
    context$param_names,
    .ctsgimme_get_cells,
    numeric(2)
  ))
  is_diagonal <- cells[, 1] == cells[, 2]
  is_group <- vapply(seq_len(nrow(cells)), function(index) {
    value <- group_drift[cells[index, 1], cells[index, 2]]
    !is.na(value) && value != "0"
  }, logical(1))
  eligible <- !is_diagonal & !is_group

  list(
    features = features[kept.ids, eligible, drop = FALSE],
    kept.ids = kept.ids,
    eligible.features = context$param_names[eligible]
  )
}

.ctsgimme_robust_scale_subgroup_features <- function(features, clip = 5) {
  center <- apply(features, 2, stats::median, na.rm = TRUE)
  spread <- apply(features, 2, stats::mad, na.rm = TRUE)
  fallback <- apply(features, 2, stats::sd, na.rm = TRUE)
  replace <- !is.finite(spread) | spread < 1e-8
  spread[replace] <- fallback[replace]
  keep <- is.finite(spread) & spread >= 1e-8
  if (!any(keep)) stop("No variable subgroup features remained after scaling.")

  scaled <- sweep(features[, keep, drop = FALSE], 2, center[keep], "-")
  scaled <- sweep(scaled, 2, spread[keep], "/")
  scaled[!is.finite(scaled)] <- 0
  scaled[scaled < -clip] <- -clip
  scaled[scaled > clip] <- clip
  list(
    scaled = scaled,
    center = center[keep],
    spread = spread[keep],
    kept.features = colnames(features)[keep],
    clip = clip
  )
}

.ctsgimme_select_recurrent_subgroup_features <- function(
    features, subject_alpha = 0.05, feature_fdr = 0.05) {
  if (!is.matrix(features)) features <- as.matrix(features)
  if (!nrow(features) || !ncol(features)) {
    stop("No subgroup score features were available for recurrence screening.")
  }
  score_threshold <- sqrt(stats::qchisq(1 - subject_alpha, df = 1))
  detected <- is.finite(features) & abs(features) >= score_threshold
  detection_count <- colSums(detected)
  recurrence_p <- stats::pbinom(
    detection_count - 1L,
    size = nrow(features),
    prob = subject_alpha,
    lower.tail = FALSE
  )
  recurrence_q <- stats::p.adjust(recurrence_p, method = "BH")
  keep <- is.finite(recurrence_q) & recurrence_q <= feature_fdr
  if (!any(keep)) {
    stop(
      "No eligible subgroup paths recurred more often than the subject-level ",
      "null rate after BH correction."
    )
  }

  signed_detection <- matrix(
    0,
    nrow = nrow(features),
    ncol = ncol(features),
    dimnames = dimnames(features)
  )
  signed_detection[detected] <- sign(features[detected])
  selected <- signed_detection[, keep, drop = FALSE]
  selection <- data.frame(
    feature = colnames(features),
    detection.count = detection_count,
    detection.prevalence = detection_count / nrow(features),
    recurrence.p = recurrence_p,
    recurrence.q = recurrence_q,
    selected = keep,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  list(
    scaled = selected,
    center = stats::setNames(rep(0, sum(keep)), colnames(selected)),
    spread = stats::setNames(rep(1, sum(keep)), colnames(selected)),
    kept.features = colnames(selected),
    clip = NA_real_,
    score.threshold = score_threshold,
    subject.alpha = subject_alpha,
    feature.fdr = feature_fdr,
    selection = selection
  )
}

.ctsgimme_subgroup_manhattan_distance <- function(scaled_features) {
  distance <- as.matrix(stats::dist(scaled_features, method = "manhattan"))
  distance / ncol(scaled_features)
}

.ctsgimme_choose_pam_by_silhouette <- function(distance_matrix, k_values) {
  n <- nrow(distance_matrix)
  k_values <- sort(unique(suppressWarnings(as.integer(k_values))))
  k_values <- k_values[is.finite(k_values) & k_values >= 2L & k_values < n]
  if (!length(k_values)) {
    stop(
      "The PAM candidate range must contain an integer in ",
      "[2, n_subjects - 1]."
    )
  }

  dissimilarity <- stats::as.dist(distance_matrix)
  fits <- lapply(k_values, function(k) {
    cluster::pam(dissimilarity, k = k, diss = TRUE)
  })
  widths <- vapply(
    fits,
    function(fit) as.numeric(fit$silinfo$avg.width),
    numeric(1)
  )
  best <- which.max(widths)
  selected <- fits[[best]]
  names(selected$clustering) <- rownames(distance_matrix)
  attr(selected, "ctsgimme.method") <- "pam"
  attr(selected, "ctsgimme.candidate.silhouettes") <- stats::setNames(
    widths,
    k_values
  )

  list(
    fit = selected,
    candidates = data.frame(
      k = k_values,
      average.silhouette = widths,
      selected = seq_along(k_values) == best,
      row.names = NULL
    )
  )
}

.ctsgimme_detect_subgroups_pam <- function(
    context, group_drift, max.subgroups) {
  evidence <- .ctsgimme_extract_subgroup_score_features(context, group_drift)
  if (length(evidence$kept.ids) < 3L) {
    return(.ctsgimme_one_subgroup_result(
      context,
      method = "pam",
      reason = paste0(
        "PAM subgrouping had fewer than three valid subject models; ",
        "assigning all subjects to one subgroup."
      ),
      warn = TRUE
    ))
  }

  scaled <- tryCatch(
    .ctsgimme_select_recurrent_subgroup_features(evidence$features),
    error = function(e) e
  )
  if (inherits(scaled, "error")) {
    return(.ctsgimme_one_subgroup_result(
      context,
      method = "pam",
      reason = paste0(
        "PAM subgrouping found no recurrent eligible features (",
        conditionMessage(scaled), "); assigning all subjects to one subgroup."
      ),
      warn = TRUE
    ))
  }

  distance <- .ctsgimme_subgroup_manhattan_distance(scaled$scaled)
  maximum_k <- min(as.integer(max.subgroups), nrow(distance) - 1L)
  selected <- .ctsgimme_choose_pam_by_silhouette(
    distance,
    seq.int(2L, maximum_k)
  )
  partial_memb <- selected$fit$clustering
  names(partial_memb) <- evidence$kept.ids
  memb <- rep(1L, length(context$ids))
  names(memb) <- as.character(context$ids)
  matched <- match(evidence$kept.ids, names(memb))
  keep <- !is.na(matched)
  memb[matched[keep]] <- as.integer(partial_memb[keep])
  if (any(!keep)) {
    warning(
      "Some PAM model IDs did not match dataframe IDs; ",
      "unmatched IDs were left in subgroup 1.",
      call. = FALSE
    )
  }

  .ctsgimme_safe_png(
    file.path(context$directory, "pam_silhouette_selection.png"),
    {
      colors <- ifelse(selected$candidates$selected, "#2166AC", "#BDBDBD")
      barplot(
        selected$candidates$average.silhouette,
        names.arg = selected$candidates$k,
        col = colors,
        ylim = c(
          min(0, selected$candidates$average.silhouette),
          max(1, selected$candidates$average.silhouette)
        ),
        xlab = "Number of subgroups (k)",
        ylab = "Average silhouette width",
        main = "PAM subgroup-count selection"
      )
    }
  )

  list(
    membership = memb,
    clustering = selected$fit,
    walktrap = NULL,
    method = "pam",
    diagnostics = list(
      raw.features = evidence$features,
      scaled.features = scaled$scaled,
      distance = distance,
      candidates = selected$candidates,
      centers = scaled$center,
      spreads = scaled$spread,
      clip = scaled$clip,
      recurrence.selection = scaled$selection,
      recurrence.score.threshold = scaled$score.threshold,
      recurrence.subject.alpha = scaled$subject.alpha,
      recurrence.feature.fdr = scaled$feature.fdr
    )
  )
}
