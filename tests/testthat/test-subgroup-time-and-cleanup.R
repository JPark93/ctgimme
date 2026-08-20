test_that("cleanup removes only shared-search intermediate files", {
  output_directory <- tempfile("ctgimme-cleanup-")
  on.exit(unlink(output_directory, recursive = TRUE, force = TRUE), add = TRUE)

  subgroup_directory <- file.path(output_directory, "Models", "Subgroup 1")
  individual_directory <- file.path(output_directory, "Models", "Individuals")
  mi_directory <- file.path(output_directory, "MIs")
  dir.create(subgroup_directory, recursive = TRUE)
  dir.create(individual_directory, recursive = TRUE)
  dir.create(mi_directory, recursive = TRUE)

  intermediate_files <- c(
    file.path(output_directory, "Models", "Model_A.RDS"),
    file.path(output_directory, "Models", "Model_B.RDS"),
    file.path(mi_directory, "MI_A.RDS"),
    file.path(mi_directory, "MI_B.RDS"),
    file.path(subgroup_directory, "Model_A.RDS"),
    file.path(subgroup_directory, "Model_B.RDS")
  )
  retained_files <- c(
    file.path(individual_directory, "FinalModel_A.RDS"),
    file.path(subgroup_directory, "Subgroup_1Model.RDS"),
    file.path(output_directory, "GStruc.RDS")
  )
  for (path in c(intermediate_files, retained_files)) {
    saveRDS(basename(path), path)
  }
  context <- list(directory = output_directory, ids = c("A", "B"))

  expect_invisible(ctgimme:::.ctgimme_cleanup(context, TRUE))
  expect_true(all(file.exists(c(intermediate_files, retained_files))))
  expect_true(dir.exists(mi_directory))

  expect_invisible(ctgimme:::.ctgimme_cleanup(context, FALSE))
  expect_false(any(file.exists(intermediate_files)))
  expect_true(all(file.exists(retained_files)))
  expect_length(list.files(mi_directory, all.files = TRUE, no.. = TRUE), 0L)
})

test_that("intermediate cleanup treats wildcard characters literally", {
  output_directory <- tempfile("ctgimme-literal-cleanup-")
  on.exit(
    unlink(output_directory, recursive = TRUE, expand = FALSE),
    add = TRUE
  )
  mi_directory <- file.path(output_directory, "MIs")
  model_directory <- file.path(output_directory, "Models")
  dir.create(mi_directory, recursive = TRUE)
  dir.create(model_directory, recursive = TRUE)

  literal_files <- c(
    file.path(mi_directory, "MI_A[1].RDS"),
    file.path(model_directory, "Model_A[1].RDS")
  )
  sibling_files <- c(
    file.path(mi_directory, "MI_A1.RDS"),
    file.path(model_directory, "Model_A1.RDS")
  )
  for (path in c(literal_files, sibling_files)) saveRDS(TRUE, path)
  context <- list(directory = output_directory)

  expect_invisible(
    ctgimme:::.ctgimme_delete_subject_files(
      context,
      "A[1]",
      model_directory
    )
  )
  expect_false(any(file.exists(literal_files)))
  expect_true(all(file.exists(sibling_files)))
})
