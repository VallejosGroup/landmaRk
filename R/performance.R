#' Performance metrics
#'
#' Computes concordance index (c-index) and Brier scores at the specified
#' landmark times and prediction horizons.
#'
#' @param x An object of class \code{\link{LandmarkAnalysis}}.
#' @param landmarks A numeric vector of landmark times.
#' @param horizons Vector of prediction horizons up to when the survival
#'   submodel is fitted.
#' @param c_index A logical. If TRUE (default), C index is reported.
#' @param brier A logical. If TRUE (default), Brier score is reported.
#' @param auc_t A logical. If TRUE (default), AUC_t is reported.
#' @param train A logical. If TRUE (default), performance metrics are computed
#'   in the training set. If FALSE, they are computed in the test set.
#'
#' @returns Data frame with performance metrics across the specified landmark
#' times and prediction horizons.
#' @export
#'
#' @examples
setGeneric(
  "performance_metrics",
  function(
    x,
    landmarks,
    horizons,
    c_index = TRUE,
    brier = TRUE,
    auc_t = FALSE,
    train = TRUE
  ) {
    standardGeneric("performance_metrics")
  }
)

#' Performance metrics
#'
#' Computes concordance index (c-index) and Brier scores at the specified
#' landmark times and prediction horizons.
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
  function(
    x,
    landmarks,
    horizons,
    c_index = TRUE,
    brier = TRUE,
    auc_t = FALSE,
    train = TRUE
  ) {
    error_str <- NULL
    fold <- NULL
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
    auct_list <- list()
    for (i in seq_along(landmarks)) {
      landmark <- landmarks[i]
      horizon <- horizons[i]
      at_risk_individuals <- x@risk_sets[[as.character(landmark)]]

      # Retrieve survival analysis dataset (censor events past horizon time)
      dataset <- x@survival_datasets[[paste0(landmark, "-", horizon)]]

      # Recover the observations and predictions (in-sample or out-of-sample)
      if (train) {
        dataset <- x@survival_datasets[[paste0(landmark, "-", horizon)]]
        predictions <- x@survival_predictions[[paste0(landmark, "-", horizon)]]
      } else {
        dataset <- x@survival_datasets_test[[paste0(landmark, "-", horizon)]]
        predictions <- x@survival_predictions_test[[paste0(
          landmark,
          "-",
          horizon
        )]]
      }

      if (brier) {
        brier_list[[paste0(landmark, "-", horizon)]] <-
          .BinaryBrierScore(
            predictions = predictions,
            time = dataset$event_time,
            status = dataset$event_status,
            tau = horizon,
            cause = 1
          )
      }
      if (c_index) {
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
      if (auc_t) {
        timepoints <- seq(min(dataset[, "event_time"]), max(dataset[, "event_time"]), length.out = 12)
        timepoints <- timepoints[-c(1,length(timepoints))]
        auct_list[[paste0(landmark, "-", horizon)]] <-
          unname(timeROC::timeROC(
            T = dataset[, "event_time"],
            delta = dataset[, "event_status"],
            marker = predictions,
            cause = 1,
            times = timepoints
          )$AUC)
      }
      # if (auc_t) {
      #   auct_list[[paste0(landmark, "-", horizon)]] <-
      #     prueba <- timeROC::timeROC(
      #       T = dataset[, "event_time"],
      #       delta = dataset[, "event_status"],
      #       marker = unname(predictions),
      #       cause = 1,
      #       weighting = "marginal",
      #       times = 365.25,
      #       iid = TRUE
      #     )
      #
      # }
    }
    if (c_index) {
      scores <- cbind(scores, cindex = unlist(cindex_list))
    }
    if (brier) {
      scores <- cbind(scores, brier = unlist(brier_list))
    }
    if (auc_t) {
      auct_matrix <- do.call(rbind, auct_list)
      colnames(auct_matrix) <- paste0("AUCt", 1:ncol(auct_matrix))
      scores <- cbind(scores, auct_matrix)
    }
    return(scores)
  }
)
