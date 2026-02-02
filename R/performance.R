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
#' @param auc_t A logical. If TRUE, AUC_t is reported.
#' @param train A logical. If TRUE (default), performance metrics are computed
#'   in the training set. If FALSE, they are computed in the test set.
#' @param h_times A numeric vector of horizon times where auc_t and Brier score
#'   are calculated.
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
    train = TRUE,
    h_times = c()
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
    train = TRUE,
    h_times = c()
  ) {
    error_str <- NULL
    model <- NULLBrier <- NULL
    Brier <- NULL
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
      if (sum(dataset$event_status) == 0) {
        stop("No events in the evaluation set")
      }

      if (brier) {
        if (length(h_times) == 0) {
          brier_list[[paste0(landmark, "-", horizon)]] <- riskRegression::Score(
            # object = list(x@survival_predictions_test[[paste0(landmark, "-", horizon)]]),
            object = list(x@survival_fits[[paste0(landmark, "-", horizon)]]),
            formula = x@survival_fits[[paste0(landmark, "-", horizon)]]$formula,
            data = dataset,
            cause = 1,
            times = horizon - landmark,
            cens.method = "ipcw",
            cens.model = "km"
          )$Brier$score |>
            filter(model != "Null model") |>
            pull(Brier)
        } else {
          brier_list[[paste0(landmark, "-", horizon)]] <- riskRegression::Score(
            # object = list(x@survival_predictions_test[[paste0(landmark, "-", horizon)]]),
            object = list(x@survival_fits[[paste0(landmark, "-", horizon)]]),
            formula = x@survival_fits[[paste0(landmark, "-", horizon)]]$formula,
            data = dataset,
            cause = 1,
            times = h_times,
            cens.method = "ipcw",
            cens.model = "km"
          )$Brier$score |>
            filter(model != "Null model") |>
            pull(Brier)

          names(brier_list[[paste0(landmark, "-", horizon)]]) <- paste0(
            "Brier(",
            landmark + h_times,
            ")"
          )
        }
      }
      if (c_index) {
        if (train == TRUE) {
          survival_dataset <- x@survival_datasets[[paste0(
            landmark,
            "-",
            horizon
          )]]
        } else {
          survival_dataset <- x@survival_datasets_test[[paste0(
            landmark,
            "-",
            horizon
          )]]
        }
        cindex_list[[paste0(landmark, "-", horizon)]] <- pec::cindex(
          list(x@survival_fits[[paste0(landmark, "-", horizon)]]),
          x@survival_fits[[paste0(landmark, "-", horizon)]]$formula,
          survival_dataset,
        )$AppCindex$coxph
      }
      if (auc_t) {
        if (length(h_times) == 0) {
          auct_list[[paste0(landmark, "-", horizon)]] <- unname(
            timeROC::timeROC(
              T = dataset[, "event_time"],
              delta = dataset[, "event_status"],
              marker = predictions,
              cause = 1,
              times = horizon - landmark
            )$AUC[2]
          )
        } else {
          auct_list[[paste0(landmark, "-", horizon)]] <- unname(
            timeROC::timeROC(
              T = dataset[, "event_time"],
              delta = dataset[, "event_status"],
              marker = predictions,
              cause = 1,
              times = h_times
            )$AUC
          )
          names(auct_list[[paste0(landmark, "-", horizon)]]) <- paste0(
            "AUC(",
            landmark + h_times,
            ")"
          )
        }
      }
    }
    if (c_index) {
      scores <- cbind(scores, cindex = unlist(cindex_list))
    }
    if (brier) {
      brier_matrix <- do.call(rbind, brier_list)
      if (ncol(brier_matrix) == 1) {
        colnames(brier_matrix) <- "Brier"
      }
      scores <- cbind(scores, brier_matrix)
    }
    if (auc_t) {
      auct_matrix <- do.call(rbind, auct_list)
      if (ncol(auct_matrix) == 1) {
        colnames(auct_matrix) <- "AUCt"
      }
      scores <- cbind(scores, auct_matrix)
    }
    return(scores)
  }
)
