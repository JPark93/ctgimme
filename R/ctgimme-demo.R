#' Run a Fast ctgimme Subgrouping Demonstration
#'
#' Demonstrates the recurrent-evidence screening, distance calculation, and
#' PAM subgroup selection used by `ctgimme()` without fitting OpenMx models.
#' The deterministic example is deliberately small and completes in well under
#' a second.
#'
#' @return A list with the synthetic signed path-score `features`, the
#'   feature-selection diagnostics, the mean Manhattan `distance` matrix,
#'   candidate silhouette widths in `candidates`, and the named subgroup
#'   `membership`.
#'
#' @examples
#' demo_result <- ctgimme_demo()
#' demo_result$membership
#' demo_result$candidates
#'
#' @seealso [ctgimme()]
#' @export
ctgimme_demo <- function() {
  features <- rbind(
    S01 = c(-3.2, -2.8, 0),
    S02 = c(-3.0, -2.7, 0),
    S03 = c(-2.9, -3.1, 0),
    S04 = c(-3.1, -2.9, 0),
    S05 = c(3.2, 2.8, 0),
    S06 = c(3.0, 2.7, 0),
    S07 = c(2.9, 3.1, 0),
    S08 = c(3.1, 2.9, 0)
  )
  colnames(features) <- c("x1_to_x2", "x2_to_x1", "null_path")

  selected <- .ctgimme_select_recurrent_subgroup_features(features)
  distance <- .ctgimme_subgroup_manhattan_distance(selected$scaled)
  pam_result <- .ctgimme_choose_pam_by_silhouette(distance, 2L)

  list(
    features = features,
    selection = selected$selection,
    distance = distance,
    candidates = pam_result$candidates,
    membership = pam_result$fit$clustering
  )
}
