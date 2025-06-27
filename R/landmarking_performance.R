#' Performance metrics
#'
#' Computes concordance index (c-index) and Brier scores at the specified landmark
#' times and prediction windows.
#'
#' @param x An object of class \code{\link{Landmarking}}.
#' @param landmarks A numeric vector of landmark times.
#' @param windows A vector of prediction windows determining horizon times.
#' @param c_index A logical. If TRUE (default), C index is reported.
#' @param brier A logical. If TRUE (default), Brier score is reported.
#'
#' @returns Data frame with performance metrics across the specified landmark
#' times and prediction windows.
#' @export
#'
#' @examples
setGeneric(
  "performance_metrics",
  function(x, landmarks, windows, c_index = TRUE, brier = TRUE) {
    standardGeneric("performance_metrics")
  }
)

#' Performance metrics
#'
#' Computes concordance index (c-index) and Brier scores at the specified landmark
#' times and prediction windows.
#'
#' @inheritParams performance_metrics
#'
#' @returns
#' @export
#'
#' @examples
setMethod(
  "performance_metrics",
  "Landmarking",
  function(x, landmarks, windows, c_index = TRUE, brier = TRUE) {
    error_str <- NULL
    if (!inherits(x, "Landmarking")) {
      error_str <- c(error_str, "@x must be an object of class Landmarking")
    }
    if (is(landmarks)[1] != "numeric") {
      error_str <- c(error_str, "@landmarks must be a vector of numeric values")
    }
    if (is(windows)[1] != "numeric") {
      error_str <- c(error_str, "@windows must be a vector of numeric values")
    }
    if (is(c_index)[1] != "logical") {
      error_str <- c(error_str, "@c_index must be a logical")
    }
    if (is(brier)[1] != "logical") {
      error_str <- c(error_str, "@brier must be a logical")
    }
    if (length(error_str) > 0) {
      stop(paste(error_str, collapse = ". "))
    }
    scores <- expand.grid(landmark = landmarks, window = windows)
    brier_list <- list()
    cindex_list <- list()
    for (landmark in landmarks) {
      for (window in windows) {
        horizon <- landmark + window
        at_risk_individuals <- x@risk_sets[[as.character(landmark)]]

        # Construct dataset for survival analysis (censor events past horizon time)
        dataset <- x@data_static[at_risk_individuals, ] |>
          mutate(
            event_status = ifelse(
              get(x@event_time) > horizon,
              0,
              get(x@event_indicator)
            ),
            event_time = ifelse(
              get(x@event_time) > horizon,
              horizon,
              get(x@event_time)
            )
          )

        predictions <- x@survival_predictions[[paste0(landmark, "-", window)]]
        if (brier == TRUE) {
          brier_list[[paste0(landmark, "-", window)]] <-
            BinaryBrierScore(
              predictions = predictions,
              time = dataset$event_time,
              status = dataset$event_status,
              tau = horizon,
              cause = 1
            )
        }
        if (c_index == TRUE) {
          cindex_list[[paste0(landmark, "-", window)]] <-
            CIndexCRisks(
              predictions = predictions,
              time = dataset$event_time,
              status = dataset$event_status,
              tau = horizon,
              cause = 1,
              method = "survival",
              cens.code = 0
            )
        }
      }
    }
    if (c_index == TRUE) {
      scores <- cbind(scores, cindex = unlist(cindex_list))
    }
    if (brier == TRUE) {
      scores <- cbind(scores, brier = unlist(brier_list))
    }
    return(scores)
  }
)


performance_metrics_cv <- function(
  x,
  landmarks,
  windows,
  K = 10,
  c_index = TRUE,
  brier = TRUE
) {
  browser()
  error_str <- NULL
  if (!inherits(x, "Landmarking")) {
    error_str <- c(error_str, "@x must be an object of class Landmarking")
  }
  if (is(landmarks)[1] != "numeric") {
    error_str <- c(error_str, "@landmarks must be a vector of numeric values")
  }
  if (is(windows)[1] != "numeric") {
    error_str <- c(error_str, "@windows must be a vector of numeric values")
  }
  if (is(c_index)[1] != "logical") {
    error_str <- c(error_str, "@c_index must be a logical")
  }
  if (is(brier)[1] != "logical") {
    error_str <- c(error_str, "@brier must be a logical")
  }
  if (length(error_str) > 0) {
    stop(paste(error_str, collapse = ". "))
  }

  x@data_static$cv_fold <- sample(K, size = nrow(x@data_static), replace = TRUE)

  num_longitudinal_fits <- length(x@longitudinal_fits)
  for (name in names(x@survival_fits)) {
    # Iterate over every combination of landmark+window
    for (k in 1:K) {
      # Iterate over crossvalidation folds
    }
    if (brier == TRUE) {
      brier_list[[paste0(landmark, "-", window)]] <-
        BinaryBrierScore(
          predictions = predictions,
          time = dataset$event_time,
          status = dataset$event_status,
          tau = horizon,
          cause = 1
        )
    }
    if (c_index == TRUE) {
      cindex_list[[paste0(landmark, "-", window)]] <-
        CIndexCRisks(
          predictions = predictions,
          time = dataset$event_time,
          status = dataset$event_status,
          tau = horizon,
          cause = 1,
          method = "survival",
          cens.code = 0
        )
    }
  }

  scores <- expand.grid(landmark = landmarks, window = windows)
  brier_list <- list()
  cindex_list <- list()
  for (landmark in landmarks) {
    for (window in windows) {
      horizon <- landmark + window
      at_risk_individuals <- x@risk_sets[[as.character(landmark)]]

      # Construct dataset for survival analysis (censor events past horizon time)
      dataset <- x@data_static[at_risk_individuals, ] |>
        mutate(
          event_status = ifelse(
            get(x@event_time) > horizon,
            0,
            get(x@event_indicator)
          ),
          event_time = ifelse(
            get(x@event_time) > horizon,
            horizon,
            get(x@event_time)
          )
        )

      predictions <- x@survival_predictions[[paste0(landmark, "-", window)]]
      if (brier == TRUE) {
        brier_list[[paste0(landmark, "-", window)]] <-
          BinaryBrierScore(
            predictions = predictions,
            time = dataset$event_time,
            status = dataset$event_status,
            tau = horizon,
            cause = 1
          )
      }
      if (c_index == TRUE) {
        cindex_list[[paste0(landmark, "-", window)]] <-
          CIndexCRisks(
            predictions = predictions,
            time = dataset$event_time,
            status = dataset$event_status,
            tau = horizon,
            cause = 1,
            method = "survival",
            cens.code = 0
          )
      }
    }
  }
  if (c_index == TRUE) {
    scores <- cbind(scores, cindex = unlist(cindex_list))
  }
  if (brier == TRUE) {
    scores <- cbind(scores, brier = unlist(brier_list))
  }
  return(scores)
}
