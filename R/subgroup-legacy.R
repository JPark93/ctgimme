.ctgimme_detect_subgroups_legacy <- function(context, sub.sig.thrsh, conduct) {
  walktrap_comm <- NULL
  if (sub.sig.thrsh == 1.00) {
    return(.ctgimme_one_subgroup_result(
      context,
      method = "legacy",
      reason = "Subgroup detection was disabled because sub.sig.thrsh = 1."
    ))
  }

  .ctgimme_inform(context$verbose, "Beginning subgrouping stage.")
  rdss1 <- mixedsort(list.files(
    file.path(context$directory, "Models"),
    pattern = "^Model_.*\\.RDS$",
    full.names = TRUE
  ))
  rdss2 <- mixedsort(list.files(
    file.path(context$directory, "MIs"),
    pattern = "^MI_.*\\.RDS$",
    full.names = TRUE
  ))
  expected_ids <- as.character(context$ids)
  model_ids <- .ctgimme_artifact_subject_id(rdss1, "Model_")
  mi_ids <- .ctgimme_artifact_subject_id(rdss2, "MI_")
  rdss1 <- rdss1[
    !is.na(model_ids) & model_ids %in% expected_ids
  ]
  rdss2 <- rdss2[
    !is.na(mi_ids) & mi_ids %in% expected_ids
  ]
  mi_ids <- .ctgimme_artifact_subject_id(rdss2, "MI_")
  mi_by_id <- setNames(rdss2, mi_ids)
  models <- list()
  feature.weights <- list()
  kept.ids <- character(0)

  common_path_weight <- 2.0
  dominant_path_weight <- 1.0
  weaker_path_weight <- 0.50
  shared_zero_weight <- 0.05
  pair.idx <- which(
    upper.tri(matrix(NA, context$nvar, context$nvar)),
    arr.ind = TRUE
  )

  for (file in seq_along(rdss1)) {
    model_id <- .ctgimme_artifact_subject_id(rdss1[file], "Model_")
    mi_file <- unname(mi_by_id[model_id])
    if (length(mi_file) == 0 || is.na(mi_file)) next
    temp1 <- tryCatch(readRDS(rdss1[file]), error = function(e) NULL)
    temp2 <- tryCatch(readRDS(mi_file), error = function(e) NULL)
    if (is.null(temp1) || is.null(temp2)) next
    drifts <- tryCatch(
      subset(summary(temp1)$parameters, matrix == "A"),
      error = function(e) NULL
    )
    if (is.null(drifts) || !nrow(drifts)) next
    if ("name" %in% names(drifts)) {
      cells <- matrix(
        as.numeric(unlist(regmatches(drifts$name, gregexpr("\\d+", drifts$name)))),
        ncol = 2,
        byrow = TRUE
      )
    } else {
      cells <- cbind(as.numeric(drifts$row), as.numeric(drifts$col))
    }
    MI_vals <- temp2$MI.Full
    EPC_vals <- temp2$EPC
    if (is.null(MI_vals) || is.null(EPC_vals) || !length(MI_vals)) next
    MI.cells <- matrix(
      as.numeric(unlist(regmatches(
        names(MI_vals),
        gregexpr("\\d+", names(MI_vals))
      ))),
      ncol = 2,
      byrow = TRUE
    )
    MI.cells <- cbind(MI.cells, MI_vals, EPC_vals)

    mi.raw.mat <- matrix(NA_real_, context$nvar, context$nvar)
    if (nrow(MI.cells) > 0) {
      for (ii in seq_len(nrow(MI.cells))) {
        r <- MI.cells[ii, 1]
        c <- MI.cells[ii, 2]
        mi.raw <- suppressWarnings(abs(as.numeric(MI.cells[ii, 3])))
        if (is.finite(r) && is.finite(c) &&
            r >= 1 && r <= context$nvar &&
            c >= 1 && c <= context$nvar && is.finite(mi.raw)) {
          mi.raw.mat[r, c] <- mi.raw
        }
      }
    }
    MI.cells[, 3] <- ifelse(
      MI.cells[, 3] > qchisq(0.975, 1),
      MI.cells[, 4],
      0
    )
    temp.mat1 <- temp.mat2 <- matrix(NA, context$nvar, context$nvar)
    temp.w1 <- matrix(common_path_weight, context$nvar, context$nvar)
    temp.w2 <- matrix(weaker_path_weight, context$nvar, context$nvar)

    for (ii in seq_len(nrow(cells))) {
      val <- drifts[ii, "Estimate"]
      se <- drifts[ii, "Std.Error"]
      if (is.finite(val) && is.finite(se) && se > 0) {
        temp.mat1[cells[ii, 1], cells[ii, 2]] <- ifelse(
          abs(val) / se > qnorm(0.95),
          val,
          0
        )
      }
    }
    for (ii in seq_len(nrow(MI.cells))) {
      temp.mat2[MI.cells[ii, 1], MI.cells[ii, 2]] <- MI.cells[ii, 3]
    }
    if (nrow(pair.idx) > 0) {
      for (pp in seq_len(nrow(pair.idx))) {
        r <- pair.idx[pp, 1]
        c <- pair.idx[pp, 2]
        mi.forward <- mi.raw.mat[r, c]
        mi.reverse <- mi.raw.mat[c, r]
        if (!is.finite(mi.forward)) mi.forward <- 0
        if (!is.finite(mi.reverse)) mi.reverse <- 0
        if (mi.forward > mi.reverse) {
          temp.w2[r, c] <- dominant_path_weight
          temp.w2[c, r] <- weaker_path_weight
        } else if (mi.reverse > mi.forward) {
          temp.w2[c, r] <- dominant_path_weight
          temp.w2[r, c] <- weaker_path_weight
        } else {
          temp.w2[r, c] <- weaker_path_weight
          temp.w2[c, r] <- weaker_path_weight
        }
      }
    }
    models[[length(models) + 1]] <- cbind(temp.mat1, temp.mat2)
    feature.weights[[length(feature.weights) + 1]] <- cbind(temp.w1, temp.w2)
    kept.ids <- c(kept.ids, model_id)
  }

  if (length(models) < 2) {
    return(.ctgimme_one_subgroup_result(
      context,
      method = "legacy",
      reason = paste0(
        "Legacy subgrouping had fewer than two valid subject models; ",
        "assigning all subjects to one subgroup."
      ),
      warn = TRUE
    ))
  }

  adj.mat <- matrix(NA, length(models), length(models))
  for (ii in seq_along(models)) {
    for (jj in seq_along(models)) {
      x <- c(models[[ii]])
      y <- c(models[[jj]])
      sx <- sign(x)
      sy <- sign(y)
      w <- (c(feature.weights[[ii]]) + c(feature.weights[[jj]])) / 2
      same.nonzero <- !is.na(sx) & !is.na(sy) & sx == sy & sx != 0
      same.zero <- !is.na(sx) & !is.na(sy) & sx == 0 & sy == 0
      nonzero.sim <- sum(w[same.nonzero], na.rm = TRUE)
      zero.sim <- shared_zero_weight * sum(same.zero, na.rm = TRUE)
      adj.mat[ii, jj] <- nonzero.sim + zero.sim
    }
  }
  jpmodmax1 <- function(x, m) {
    m2 <- m
    threshold_value <- quantile(
      m2[upper.tri(m2, diag = FALSE)],
      x,
      na.rm = TRUE
    )
    m2[m2 <= threshold_value] <- 0
    g <- igraph::graph_from_adjacency_matrix(
      m2,
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    )
    p <- igraph::cluster_walktrap(g)
    mem <- igraph::membership(p)
    n <- length(mem)
    sk <- numeric(n)
    ask <- numeric(n)
    for (i in seq_len(n)) {
      internal_i <- 0
      external_i <- 0
      for (j in seq_len(n)) {
        if (m2[i, j] != 0) {
          if (mem[i] == mem[j]) internal_i <- internal_i + m2[i, j]
          else external_i <- external_i + m2[i, j]
        }
      }
      sk[i] <- external_i
      ask[i] <- internal_i + external_i
    }
    temp.1 <- cbind(mem = mem, sk = sk, ask = ask)
    cluster_ids <- sort(unique(mem))
    c_vals <- numeric(length(cluster_ids))
    for (k in seq_along(cluster_ids)) {
      cl <- cluster_ids[k]
      cluster_nodes <- which(mem == cl)
      cut_edges <- sum(temp.1[cluster_nodes, "sk"])
      total_edges <- sum(temp.1[cluster_nodes, "ask"])
      cluster_size <- length(cluster_nodes)
      if (total_edges == 0) c_vals[k] <- 1
      else c_vals[k] <- cut_edges / total_edges
      if (cluster_size == 1) c_vals[k] <- max(c_vals[k], 1)
    }
    conductance <- 1 - mean(c_vals)
    -1 * conductance
  }

  if (isTRUE(conduct)) {
    res <- nloptr::nloptr(
      x0 = 0.5,
      m = adj.mat,
      eval_f = jpmodmax1,
      lb = .01,
      ub = .99,
      opts = list(
        "algorithm" = "NLOPT_GN_DIRECT_L",
        "ftol_rel" = 1.0e-8,
        maxeval = 1000
      )
    )
    adj.mat[which(adj.mat <= quantile(
      adj.mat[upper.tri(adj.mat, diag = FALSE)],
      res$solution
    ))] <- 0
  } else {
    adj.mat <- as.matrix(adj.mat)
    adj.mat <- adj.mat - min(adj.mat, na.rm = TRUE)
  }
  diag(adj.mat) <- 0
  g <- igraph::graph_from_adjacency_matrix(
    adj.mat,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  walktrap_comm <- igraph::cluster_walktrap(g)
  .ctgimme_safe_png(
    file.path(context$directory, "walktrap_community_plot.png"),
    {
      plot(
        walktrap_comm,
        g,
        vertex.size = 15,
        vertex.label = kept.ids,
        main = "Walktrap Community Detection"
      )
    }
  )
  partial_memb <- igraph::membership(walktrap_comm)
  names(partial_memb) <- kept.ids
  memb <- rep(1, length(context$ids))
  names(memb) <- as.character(context$ids)
  matched <- match(kept.ids, names(memb))
  keep <- !is.na(matched)
  memb[matched[keep]] <- partial_memb[keep]
  if (any(!keep)) {
    warning(
      "Some subgroup model IDs did not match dataframe IDs; ",
      "unmatched IDs were left in subgroup 1."
    )
  }

  attr(walktrap_comm, "ctgimme.method") <- "legacy"
  list(
    membership = memb,
    clustering = walktrap_comm,
    walktrap = walktrap_comm,
    method = "legacy",
    diagnostics = list(
      similarity = adj.mat,
      kept.ids = kept.ids,
      conductance.optimized = isTRUE(conduct)
    )
  )
}
