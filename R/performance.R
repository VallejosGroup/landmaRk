#' Performance metrics
#'
#' Computes concordance index (c-index) and Brier scores at the specified landmark
#' times and prediction horizons.
#'
#' @param x An object of class \code{\link{LandmarkAnalysis}}.
#' @param landmarks A numeric vector of landmark times.
#' @param horizons Vector of prediction horizons up to when the survival submodel
#'   is fitted.
#' @param c_index A logical. If TRUE (default), C index is reported.
#' @param brier A logical. If TRUE (default), Brier score is reported.
#'
#' @returns Data frame with performance metrics across the specified landmark
#' times and prediction horizons.
#' @export
#'
#' @examples
setGeneric(
  "performance_metrics",
  function(x, landmarks, horizons, c_index = TRUE, brier = TRUE) {
    standardGeneric("performance_metrics")
  }
)

#' Performance metrics
#'
#' Computes concordance index (c-index) and Brier scores at the specified landmark
#' times and prediction horizons.
#'
#' @inheritParams performance_metrics
#'
#' @returns
#' @export
#'
#' @examples
setMethod(
  "performance_metrics",
  "LandmarkAnalysis",
  function(x, landmarks, horizons, c_index = TRUE, brier = TRUE) {
    error_str <- NULL
    if (!inherits(x, "LandmarkAnalysis")) {
      error_str <- c(
        error_str,
        "@x must be an object of class LandmarkAnalysis"
      )
    }
    if (is(landmarks)[1] != "numeric") {
      error_str <- c(error_str, "@landmarks must be a vector of numeric values")
    }
    if (is(horizons)[1] != "numeric") {
      error_str <- c(error_str, "@horizons must be a vector of numeric values")
    }
    if (is(c_index)[1] != "logical") {
      error_str <- c(error_str, "@c_index must be a logical")
    }
    if (is(brier)[1] != "logical") {
      error_str <- c(error_str, "@brier must be a logical")
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

    scores <- cbind(landmark = landmarks, horizon = horizons)
    brier_list <- list()
    cindex_list <- list()
    for (i in seq_along(landmarks)) {
      landmark <- landmarks[i]
      horizon <- horizons[i]
      at_risk_individuals <- x@risk_sets[[as.character(landmark)]]

      # Construct dataset for survival analysis (censor events past horizon time)
      dataset <- data.frame(at_risk_individuals)
      colnames(dataset) <- x@ids
      dataset <- dataset |>
        left_join(x@data_static, by = stats::setNames(x@ids, x@ids)) |>
        mutate(
          event_status = ifelse(
            get(x@event_time) > horizon,
            0,
            get(x@event_indicator)
          ),
          event_time = ifelse(
            get(x@event_time) > horizon,
            horizon - landmark,
            get(x@event_time) - landmark
          )
        )

      predictions <- x@survival_predictions[[paste0(landmark, "-", horizon)]]
      if (brier == TRUE) {
        brier_list[[paste0(landmark, "-", horizon)]] <-
          .BinaryBrierScore(
            predictions = predictions,
            time = dataset$event_time,
            status = dataset$event_status,
            tau = horizon,
            cause = 1
          )
      }
      if (c_index == TRUE) {
        cindex_list[[paste0(landmark, "-", horizon)]] <-
          .CIndexCRisks(
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
    if (c_index == TRUE) {
      scores <- cbind(scores, cindex = unlist(cindex_list))
    }
    if (brier == TRUE) {
      scores <- cbind(scores, brier = unlist(brier_list))
    }
    return(scores)
  }
)
