.ctgimme_cleanup <- function(context, keep.intermediate) {
  if (!keep.intermediate) {
    .ctgimme_delete_subject_files(
      context,
      context$ids,
      file.path(context$directory, "Models")
    )
    subgroup_dirs <- list.dirs(
      file.path(context$directory, "Models"),
      recursive = FALSE,
      full.names = TRUE
    )
    subgroup_dirs <- subgroup_dirs[
      grepl("Subgroup ", basename(subgroup_dirs), fixed = TRUE)
    ]
    for (sd in subgroup_dirs) {
      unlink(
        file.path(sd, paste0("Model_", as.character(context$ids), ".RDS")),
        force = TRUE,
        expand = FALSE
      )
    }

    mi_dir <- file.path(context$directory, "MIs")
    if (dir.exists(mi_dir) &&
        !length(list.files(mi_dir, all.files = TRUE, no.. = TRUE))) {
      unlink(mi_dir, recursive = FALSE, force = TRUE, expand = FALSE)
    }
  }
  invisible(NULL)
}
