#' Continuous-Time Subgrouping with GIMME
#'
#' `ctgimme` estimates group-, subgroup-, and individual-level dynamic
#' structures using continuous-time state-space models. The package uses a
#' shrunk empirical initial-state covariance and recurrent-evidence PAM as its
#' default subgrouping procedure.
#'
#' @references Park, J. J., Fisher, Z. F., Hunter, M. D., Shenk, C., Russell,
#'   M., Molenaar, P. C. M., & Chow, S.-M. (2025). Unsupervised model
#'   construction in continuous-time. *Structural Equation Modeling: A
#'   Multidisciplinary Journal, 32*(3), 377--399.
#'   \doi{10.1080/10705511.2024.2429544}
#'
#' @import OpenMx
#' @importFrom graphics barplot plot
#' @importFrom grDevices dev.off png
#' @importFrom gtools mixedsort
#' @importFrom methods is
#' @importFrom parallel clusterEvalQ clusterExport makeCluster parLapply
#'   stopCluster
#' @importFrom qgraph qgraph
#' @importFrom stats as.dist coef cov dist mad median p.adjust pbinom qchisq
#'   pchisq qnorm quantile sd setNames
#' @keywords internal
"_PACKAGE"
