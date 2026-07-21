# Validate inputs for performance_metrics function
.validate_performance_inputs <- function(
  x,
  landmarks,
  horizons,
  c_index,
  brier,
  auc_t,
  train,
  h_times,
  cause
) {
  error_str <- NULL

  if (!inherits(x, "LandmarkAnalysis")) {
    error_str <- c(
      error_str,
      "@x must be an object of class LandmarkAnalysis"
    )
  }
  if (!is.numeric(landmarks)) {
    error_str <- c(error_str, "@landmarks must be a vector of numeric values")
  }
  if (!is.numeric(horizons)) {
    error_str <- c(error_str, "@horizons must be a vector of numeric values")
  }
  if (!(is.logical(c_index) && length(c_index) == 1L && !is.na(c_index))) {
    error_str <- c(error_str, "@c_index must be a length-1, non-NA logical")
  }
  if (!(is.logical(brier) && length(brier) == 1L && !is.na(brier))) {
    error_str <- c(error_str, "@brier must be a length-1, non-NA logical")
  }
  if (!(is.logical(auc_t) && length(auc_t) == 1L && !is.na(auc_t))) {
    error_str <- c(error_str, "@auc_t must be a length-1, non-NA logical")
  }
  if (!(is.logical(train) && length(train) == 1L && !is.na(train))) {
    error_str <- c(error_str, "@train must be a length-1, non-NA logical")
  }
  if (
    !(is.numeric(cause) &&
      length(cause) == 1L &&
      !is.na(cause) &&
      is.finite(cause) &&
      cause > 0 &&
      cause %% 1 == 0)
  ) {
    error_str <- c(
      error_str,
      "@cause must be a length-1, finite, positive integer"
    )
  }
  if (length(landmarks) != length(horizons)) {
    error_str <- c(
      error_str,
      "@landmarks and @horizons must be of the same length"
    )
  }

  if (length(error_str) > 0) {
    stop(paste(error_str, collapse = ". "))
  }
}
