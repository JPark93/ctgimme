test_that("recurrent evidence retains repeated signed detections", {
  features <- cbind(
    recurrent = rep(c(-3, 3), 10),
    null = rep(0, 20)
  )
  selected <- ctsgimme:::.ctsgimme_select_recurrent_subgroup_features(features)

  expect_identical(selected$kept.features, "recurrent")
  expect_identical(as.numeric(selected$scaled[, 1]), rep(c(-1, 1), 10))
  expect_true(selected$selection$selected[1])
  expect_false(selected$selection$selected[2])
})

test_that("subgroup distance is mean Manhattan distance per feature", {
  features <- rbind(a = c(-1, 0), b = c(1, 1), c = c(-1, 1))
  distance <- ctsgimme:::.ctsgimme_subgroup_manhattan_distance(features)

  expect_equal(distance["a", "b"], 1.5)
  expect_equal(distance["a", "c"], 0.5)
  expect_equal(distance["b", "c"], 1)
})

test_that("PAM selection preserves subject identifiers and diagnostics", {
  points <- rbind(
    s1 = c(0, 0), s2 = c(0, 0.1),
    s3 = c(4, 4), s4 = c(4, 4.1),
    s5 = c(8, 0), s6 = c(8, 0.1)
  )
  distance <- as.matrix(stats::dist(points, method = "manhattan"))
  selected <- ctsgimme:::.ctsgimme_choose_pam_by_silhouette(distance, 2:3)

  expect_identical(names(selected$fit$clustering), rownames(points))
  expect_identical(attr(selected$fit, "ctsgimme.method"), "pam")
  expect_equal(nrow(selected$candidates), 2L)
  expect_equal(sum(selected$candidates$selected), 1L)
})

test_that("subgroup membership artifacts identify every subject", {
  output_directory <- tempfile("ctsgimme-membership-")
  membership <- stats::setNames(c(2L, 1L, 2L, 1L), c("S04", "S01", "S08", "S03"))
  features <- rbind(
    S04 = c(4, 4),
    S01 = c(0, 0),
    S08 = c(4.2, 4.1),
    S03 = c(0.2, 0.1)
  )
  distance <- as.matrix(stats::dist(features, method = "manhattan"))

  artifacts <- ctsgimme:::.ctsgimme_write_subgroup_membership(
    membership,
    output_directory,
    distance = distance,
    method = "pam"
  )
  written <- utils::read.csv(artifacts$membership, stringsAsFactors = FALSE)

  expect_true(file.exists(artifacts$plot))
  expect_gt(file.info(artifacts$plot)$size, 0)
  expect_identical(basename(artifacts$plot), "Subgroups Plot.png")
  expect_identical(artifacts$plot.type, "distance-map")
  expect_identical(artifacts$distance.source, "PAM mean Manhattan distance")
  expect_identical(written$subject, names(membership))
  expect_identical(written$subgroup, unname(membership))
})

test_that("subject map coordinates preserve pairwise proximity", {
  points <- rbind(
    A = c(0, 0),
    B = c(1, 0),
    C = c(0, 2),
    D = c(3, 2)
  )
  distance <- as.matrix(stats::dist(points))
  map <- ctsgimme:::.ctsgimme_subject_distance_coordinates(distance)
  mapped_distance <- as.matrix(stats::dist(map$points))

  expect_gt(stats::cor(distance[upper.tri(distance)],
                       mapped_distance[upper.tri(mapped_distance)]), 0.99)
  expect_gt(map$goodness.of.fit, 0.99)
})

test_that("membership artifacts fall back to a roster without distances", {
  output_directory <- tempfile("ctsgimme-membership-roster-")
  membership <- stats::setNames(c(1L, 1L, 2L), c("A", "B", "C"))

  artifacts <- ctsgimme:::.ctsgimme_write_subgroup_membership(
    membership,
    output_directory
  )

  expect_identical(artifacts$plot.type, "membership-roster")
  expect_null(artifacts$distance.source)
  expect_true(file.exists(artifacts$plot))
})
