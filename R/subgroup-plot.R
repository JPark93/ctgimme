# ============================================================================
# Subgroup membership artifacts and distance-map visualization
# ============================================================================

.ctsgimme_subject_distance_coordinates <- function(distance) {
  distance <- as.matrix(distance)
  subject_count <- nrow(distance)
  subject_ids <- rownames(distance)
  if (is.null(subject_ids)) subject_ids <- as.character(seq_len(subject_count))

  if (subject_count == 1L) {
    points <- matrix(c(0, 0), nrow = 1L, dimnames = list(subject_ids, NULL))
    return(list(points = points, goodness.of.fit = NA_real_, tied = TRUE))
  }

  maximum_distance <- max(distance)
  if (maximum_distance <= sqrt(.Machine$double.eps)) {
    angles <- seq(0, 2 * pi, length.out = subject_count + 1L)[-1L]
    points <- cbind(cos(angles), sin(angles)) * 0.05
    rownames(points) <- subject_ids
    return(list(points = points, goodness.of.fit = NA_real_, tied = TRUE))
  }

  dimensions <- min(2L, subject_count - 1L)
  fit <- suppressWarnings(stats::cmdscale(
    stats::as.dist(distance),
    k = dimensions,
    eig = TRUE,
    add = TRUE
  ))
  points <- as.matrix(fit$points)
  if (ncol(points) == 1L) points <- cbind(points, 0)
  rownames(points) <- subject_ids

  plot_span <- max(
    diff(range(points[, 1L])),
    diff(range(points[, 2L])),
    sqrt(.Machine$double.eps)
  )
  position_keys <- paste(
    round(points[, 1L] / plot_span, 8L),
    round(points[, 2L] / plot_span, 8L),
    sep = ":"
  )
  tied_positions <- split(seq_len(subject_count), position_keys)
  tied_positions <- tied_positions[lengths(tied_positions) > 1L]
  if (length(tied_positions)) {
    jitter_radius <- 0.025 * plot_span
    for (indices in tied_positions) {
      angles <- seq(0, 2 * pi, length.out = length(indices) + 1L)[-1L]
      points[indices, 1L] <- points[indices, 1L] + jitter_radius * cos(angles)
      points[indices, 2L] <- points[indices, 2L] + jitter_radius * sin(angles)
    }
  }

  goodness_of_fit <- if (length(fit$GOF)) {
    unname(utils::tail(fit$GOF, 1L))
  } else {
    NA_real_
  }
  list(
    points = points,
    goodness.of.fit = goodness_of_fit,
    tied = length(tied_positions) > 0L
  )
}

.ctsgimme_prepare_subject_distance <- function(
    membership, distance = NULL, similarity = NULL, diagnostic_ids = NULL) {
  source <- NULL
  subject_distance <- NULL

  if (!is.null(distance)) {
    subject_distance <- as.matrix(distance)
    source <- "PAM mean Manhattan distance"
  } else if (!is.null(similarity)) {
    similarity <- as.matrix(similarity)
    if (length(similarity) && all(is.finite(similarity))) {
      similarity <- (similarity + t(similarity)) / 2
      subject_distance <- max(similarity) - similarity
      diag(subject_distance) <- 0
      source <- "legacy model-similarity distance"
    }
  }

  if (is.null(subject_distance) || !is.numeric(subject_distance) ||
      length(dim(subject_distance)) != 2L ||
      nrow(subject_distance) != ncol(subject_distance) ||
      nrow(subject_distance) < 2L || any(!is.finite(subject_distance))) {
    return(NULL)
  }

  subject_distance <- (subject_distance + t(subject_distance)) / 2
  tolerance <- sqrt(.Machine$double.eps)
  subject_distance[subject_distance < 0 & subject_distance > -tolerance] <- 0
  if (any(subject_distance < 0)) return(NULL)
  diag(subject_distance) <- 0

  distance_ids <- rownames(subject_distance)
  column_ids <- colnames(subject_distance)
  if (is.null(distance_ids) && length(diagnostic_ids) == nrow(subject_distance)) {
    distance_ids <- as.character(diagnostic_ids)
  }
  if (is.null(column_ids) && !is.null(distance_ids)) column_ids <- distance_ids
  if (is.null(distance_ids) && !is.null(column_ids)) distance_ids <- column_ids
  if (is.null(distance_ids) && nrow(subject_distance) == length(membership)) {
    distance_ids <- names(membership)
    column_ids <- distance_ids
  }
  if (is.null(distance_ids) || is.null(column_ids) ||
      anyDuplicated(distance_ids) || anyDuplicated(column_ids) ||
      !setequal(distance_ids, column_ids)) {
    return(NULL)
  }

  rownames(subject_distance) <- distance_ids
  colnames(subject_distance) <- column_ids
  subject_distance <- subject_distance[distance_ids, distance_ids, drop = FALSE]
  mapped_ids <- intersect(names(membership), distance_ids)
  if (length(mapped_ids) < 2L) return(NULL)

  list(
    distance = subject_distance[mapped_ids, mapped_ids, drop = FALSE],
    source = source,
    omitted = setdiff(names(membership), mapped_ids)
  )
}

.ctsgimme_write_subgroup_membership <- function(
    membership, directory, distance = NULL, similarity = NULL,
    diagnostic_ids = NULL, method = NULL) {
  if (!length(membership)) {
    stop("Cannot plot subgroup membership because no assignments were supplied.")
  }
  if (anyNA(membership)) {
    stop("Cannot plot subgroup membership because some assignments are missing.")
  }

  subject_ids <- names(membership)
  if (is.null(subject_ids) || length(subject_ids) != length(membership) ||
      anyNA(subject_ids) || any(!nzchar(subject_ids))) {
    subject_ids <- as.character(seq_along(membership))
    names(membership) <- subject_ids
  }
  membership_data <- data.frame(
    subject = as.character(subject_ids),
    subgroup = unname(membership),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  csv_file <- file.path(directory, "subgroup_membership.csv")
  utils::write.csv(membership_data, csv_file, row.names = FALSE)

  subgroup_levels <- sort(unique(membership_data$subgroup))
  subgroup_ids <- lapply(subgroup_levels, function(value) {
    membership_data$subject[membership_data$subgroup == value]
  })
  subgroup_sizes <- lengths(subgroup_ids)
  n_subgroups <- length(subgroup_levels)
  subgroup_colors <- grDevices::hcl.colors(n_subgroups, palette = "Dark 3")
  names(subgroup_colors) <- as.character(subgroup_levels)

  prepared_distance <- .ctsgimme_prepare_subject_distance(
    membership,
    distance = distance,
    similarity = similarity,
    diagnostic_ids = diagnostic_ids
  )
  plot_file <- file.path(directory, "Subgroups Plot.png")

  if (!is.null(prepared_distance)) {
    mapped_ids <- rownames(prepared_distance$distance)
    mapped_membership <- membership[mapped_ids]
    map <- .ctsgimme_subject_distance_coordinates(prepared_distance$distance)
    coordinates <- map$points
    plot_span <- max(
      diff(range(coordinates[, 1L])),
      diff(range(coordinates[, 2L])),
      1
    )
    padding <- 0.16 * plot_span
    x_limits <- range(coordinates[, 1L]) + c(-padding, padding)
    y_limits <- range(coordinates[, 2L]) + c(-padding, padding)

    subject_count <- nrow(coordinates)
    neighbors <- min(3L, subject_count - 1L)
    edge_pairs <- do.call(rbind, lapply(seq_len(subject_count), function(index) {
      ordered <- order(prepared_distance$distance[index, ], na.last = NA)
      ordered <- ordered[ordered != index]
      if (!length(ordered)) return(NULL)
      cbind(index, utils::head(ordered, neighbors))
    }))
    edge_pairs <- t(apply(edge_pairs, 1L, sort))
    edge_pairs <- edge_pairs[!duplicated(
      paste(edge_pairs[, 1L], edge_pairs[, 2L], sep = ":")
    ), , drop = FALSE]
    edge_distances <- prepared_distance$distance[edge_pairs]
    maximum_distance <- max(prepared_distance$distance)
    edge_closeness <- if (maximum_distance > 0) {
      1 - edge_distances / maximum_distance
    } else {
      rep(1, length(edge_distances))
    }

    .ctsgimme_safe_png(
      plot_file,
      {
        previous_parameters <- graphics::par(
          mar = c(4.5, 1.5, 4, 1.5),
          xpd = NA
        )
        on.exit(graphics::par(previous_parameters), add = TRUE)
        graphics::plot.new()
        graphics::plot.window(
          xlim = x_limits,
          ylim = y_limits,
          asp = 1,
          xaxs = "i",
          yaxs = "i"
        )

        hull_radius <- 0.055 * plot_span
        hull_angles <- seq(0, 2 * pi, length.out = 48L)
        for (subgroup in subgroup_levels) {
          indices <- which(mapped_membership == subgroup)
          if (!length(indices)) next
          cloud <- do.call(rbind, lapply(indices, function(index) {
            cbind(
              coordinates[index, 1L] + hull_radius * cos(hull_angles),
              coordinates[index, 2L] + hull_radius * sin(hull_angles)
            )
          }))
          hull <- grDevices::chull(cloud)
          color <- subgroup_colors[[as.character(subgroup)]]
          graphics::polygon(
            cloud[hull, , drop = FALSE],
            col = grDevices::adjustcolor(color, alpha.f = 0.16),
            border = color,
            lwd = 2
          )
        }

        for (edge_index in seq_len(nrow(edge_pairs))) {
          from <- edge_pairs[edge_index, 1L]
          to <- edge_pairs[edge_index, 2L]
          closeness <- edge_closeness[edge_index]
          graphics::segments(
            coordinates[from, 1L],
            coordinates[from, 2L],
            coordinates[to, 1L],
            coordinates[to, 2L],
            col = grDevices::adjustcolor(
              "#334155",
              alpha.f = 0.15 + 0.5 * closeness
            ),
            lwd = 0.7 + 2 * closeness
          )
        }

        node_colors <- subgroup_colors[as.character(mapped_membership)]
        graphics::points(
          coordinates[, 1L],
          coordinates[, 2L],
          pch = 21,
          bg = node_colors,
          col = "#172033",
          lwd = 1.35,
          cex = 3.7
        )
        label_cex <- if (max(nchar(mapped_ids)) > 4L) 0.7 else 0.9
        graphics::text(
          coordinates[, 1L],
          coordinates[, 2L],
          labels = mapped_ids,
          cex = label_cex,
          font = 2
        )

        method_label <- if (identical(method, "pam")) "PAM" else "Subgroup"
        graphics::title(
          main = paste(method_label, "Subject Distance Map"),
          sub = "Closer subjects have more similar individual model evidence",
          line = 1
        )
        graphics::legend(
          "topright",
          legend = paste0("Subgroup ", subgroup_levels, " (n = ",
                          subgroup_sizes, ")"),
          pt.bg = subgroup_colors,
          col = "#172033",
          pch = 21,
          pt.cex = 1.7,
          bty = "n",
          cex = 0.85,
          inset = 0.01
        )

        fit_label <- if (is.finite(map$goodness.of.fit)) {
          paste0("; 2D fit = ", format(round(map$goodness.of.fit, 3), nsmall = 3))
        } else {
          ""
        }
        graphics::mtext(
          paste0(prepared_distance$source, fit_label,
                 ". Lines connect each subject to its nearest neighbors."),
          side = 1,
          line = 2.7,
          cex = 0.75,
          col = "#475569"
        )
        if (length(prepared_distance$omitted)) {
          graphics::mtext(
            paste0("No diagnostic distance for: ",
                   paste(prepared_distance$omitted, collapse = ", ")),
            side = 3,
            line = 0.15,
            cex = 0.68,
            col = "#9A3412"
          )
        }
      },
      width = 1000L,
      height = 850L
    )
    plot_type <- "distance-map"
    distance_source <- prepared_distance$source
  } else {
    max_size <- max(subgroup_sizes)
    plot_width <- max(800L, 260L * n_subgroups)
    plot_height <- max(700L, 36L * max_size + 180L)
    .ctsgimme_safe_png(
      plot_file,
      {
        graphics::plot.new()
        graphics::plot.window(
          xlim = c(0.5, n_subgroups + 0.5),
          ylim = c(0, max_size + 3),
          xaxs = "i",
          yaxs = "i"
        )
        graphics::title(main = "Subgroup Memberships", line = 1)
        for (index in seq_along(subgroup_levels)) {
          color <- subgroup_colors[index]
          graphics::rect(
            index - 0.43, 0.35, index + 0.43, max_size + 2.2,
            col = grDevices::adjustcolor(color, alpha.f = 0.14),
            border = color,
            lwd = 2
          )
          graphics::text(
            index, max_size + 1.65,
            labels = paste("Subgroup", subgroup_levels[index]),
            font = 2,
            cex = 1.05
          )
          graphics::text(
            index, max_size + 0.95,
            labels = paste0("n = ", subgroup_sizes[index]),
            cex = 0.9
          )
          ids <- subgroup_ids[[index]]
          y <- rev(seq_along(ids))
          graphics::points(
            rep(index - 0.27, length(ids)), y,
            pch = 21, bg = color, col = color, cex = 1.15
          )
          graphics::text(
            rep(index - 0.19, length(ids)), y,
            labels = ids, adj = c(0, 0.5), cex = 0.9
          )
        }
      },
      width = plot_width,
      height = plot_height
    )
    plot_type <- "membership-roster"
    distance_source <- NULL
  }

  invisible(list(
    plot = plot_file,
    membership = csv_file,
    plot.type = plot_type,
    distance.source = distance_source
  ))
}
