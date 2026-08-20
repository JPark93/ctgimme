# Aggregate modification indices
aggregate_mis <- function(rdss, param_names) {
  if (length(rdss) == 0) {
    return(list(files = NULL, EPCs = NULL))
  }
  n = length(param_names)
  nfiles = length(rdss)
  files = matrix(NA_real_, nrow = n, ncol = nfiles)
  rownames(files) = param_names
  EPCs = matrix(NA_real_, nrow = n, ncol = nfiles)
  rownames(EPCs) = param_names
  for (idx in seq_along(rdss)) {
    file = rdss[idx]
    mi_full = abs(safe_read_vector(file, "MI.Full", param_names))
    epc = safe_read_vector(file, "EPC", param_names)
    files[, idx] = mi_full
    EPCs[, idx] = epc
  }
  return(list(files = files, EPCs = EPCs))
}
