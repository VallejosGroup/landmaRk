#' Fits the specified survival model at the landmark times and up to the horizon
#' times specified by the user
#'
#' @details
#'  ## Mathematical formulation
#'  This function estimates the conditional probability of survival to horizon
#'  \eqn{s+w}, conditioned on having survived to the landmark time, \eqn{s},
#'  that is
#'  \deqn{\pi_i(s+w \vert s) = P(T_i > s+w \vert T_i \ge s, \bar{x}_i(s)), }
#'  where \eqn{i} denotes an individual's index, \eqn{T_i} is the time to event
#'  outcome for individual \eqn{i} and \eqn{\bar{x}_i(s)} are the covariates
#'  observed for individual \eqn{i}, including the observed history of dynamic
#'  covariates.
#'
#' @param x An object of class \code{\link{LandmarkAnalysis}}.
#' @param landmarks Numeric vector of landmark times.
#' @param formula A formula to be used in survival sub-model fitting.
#' @param horizons Vector of prediction horizons up to when the survival
#'   submodel is fitted.
#' @param method Method for survival analysis, either "survfit" or "coxph".
#' @param dynamic_covariates Vector of time-varying covariates to be used
#'   in the survival model.
#' @param include_clusters Boolean indicating whether to propagate cluster
#'   membership to survival analysis.
#' @param validation_fold If positive, cross-validation fold where model is
#'   fitted. If 0 (default), model fitting is performed on the complete dataset.
#'
#' @returns An object of class \code{\link{LandmarkAnalysis}}.
#' @export
#'
#' @examples
setGeneric(
  "fit_survival",
  function(
    x,
    formula,
    landmarks,
    horizons,
    method,
    dynamic_covariates = c(),
    include_clusters = FALSE,
    validation_fold = 0
  ) {
    standardGeneric("fit_survival")
  }
)

#' Fits the specified survival model at the landmark times and up to the horizon
#' times specified by the user
#'
#' @details
#'  ## Mathematical formulation
#'  This function estimates the conditional probability of survival to horizon
#'  \eqn{s+w}, conditioned on having survived to the landmark time, \eqn{s},
#'  that is
#'  \deqn{\pi_i(s+w \vert s) = P(T_i > s+w \vert T_i \ge s, \bar{x}_i(s)), }
#'  where \eqn{i} denotes an individual's index, \eqn{T_i} is the time to event
#'  outcome for individual \eqn{i} and \eqn{\bar{x}_i(s)} are the covariates
#'  observed for individual \eqn{i}, including the observed history of dynamic
#'  covariates.
#'
#' @inheritParams fit_survival
#'
#' @returns An object of class \code{\link{LandmarkAnalysis}}.
#' @export
#'
#' @examples
setMethod(
  "fit_survival",
  "LandmarkAnalysis",
  function(
    x,
    formula,
    landmarks,
    horizons,
    method,
    dynamic_covariates = c(),
    include_clusters = FALSE,
    validation_fold = 0
  ) {
    # Check that method is a function with arguments formula, data, ...
    method <- .check_method_survival_predict(method)

    if (length(landmarks) != length(horizons)) {
      stop("@landmarks and @horizons must be of the same length.")
    } else if (length(landmarks) == 1) {
      # Base case for recursion
      .check_riskset_survival(x, landmarks)

      # Recover risk sets (ids of individuals who are at risk at landmark time)
      # Construct dataset for survival analysis (censor events past horizon time)
      x@survival_datasets[[paste0(landmarks, "-", horizons)]] <-
        .create_survival_dataframe(
          x,
          landmarks,
          horizons,
          dynamic_covariates,
          include_clusters,
          validation_fold,
          train = TRUE
        )

      # Include predicted cluster membership in the training dataset and in the survival formula
      if (length(dynamic_covariates) > 0 && include_clusters == TRUE) {
        for (dynamic_covariate in dynamic_covariates) {
          formula <- as.formula(
            paste(
              as.character(formula)[2],
              as.character(formula)[1],
              paste0(
                as.character(formula)[3],
                " + ",
                paste0("cluster_", dynamic_covariate)
              )
            )
          )
        }
      }

      # Call to method that performs survival analysis
      if (is(method)[1] == "character" && method == "coxph") {
        x@survival_fits[[paste0(landmarks, "-", horizons)]] <- survival::coxph(
          formula,
          data = x@survival_datasets[[paste0(landmarks, "-", horizons)]],
          x = TRUE,
          model = TRUE
        )
      } else {
        x@survival_fits[[paste0(landmarks, "-", horizons)]] <- method(
          formula,
          data = x@survival_datasets[[paste0(landmarks, "-", horizons)]]
        )
      }
    } else {
      # Recursion
      x <- fit_survival(
        x,
        formula,
        landmarks[1],
        horizons[1],
        method,
        dynamic_covariates,
        include_clusters,
        validation_fold
      )
      x <- fit_survival(
        x,
        formula,
        landmarks[-1],
        horizons[-1],
        method,
        dynamic_covariates,
        include_clusters,
        validation_fold
      )
    }
    x
  }
)

#' Make predictions for time-to-event outcomes at specified horizon times
#'
#' @param x An object of class \code{\link{LandmarkAnalysis}}.
#' @param landmarks A numeric vector of landmark times.
#' @param horizons Vector of prediction horizons up to when the survival submodel
#'   is fitted.
#' @param method R function that is used to make predictions
#' @param dynamic_covariates Vector of time-varying covariates to be used
#'   in the survival model.
#' @param include_clusters Boolean indicating whether to propagate cluster
#'   membership to survival analysis.
#' @param validation_fold If positive, cross-validation fold where model is
#'   fitted. If 0 (default), model fitting is performed on the complete dataset.
#' @param ... Additional arguments passed to the prediction function (e.g.
#'   number of classes/clusters for lcmm).
#'
#' @returns An object of class \code{\link{LandmarkAnalysis}}.
#' @export
#'
#' @examples
setGeneric(
  "predict_survival",
  function(
    x,
    landmarks,
    horizons,
    method,
    dynamic_covariates = c(),
    include_clusters = FALSE,
    validation_fold = 0,
    ...
  ) {
    standardGeneric("predict_survival")
  }
)

#' Make predictions for time-to-event outcomes at specified horizon times
#'
#' @inheritParams predict_survival
#'
#' @returns An object of class \code{\link{LandmarkAnalysis}}.
#' @export
#'
#' @examples
setMethod(
  "predict_survival",
  "LandmarkAnalysis",
  function(
    x,
    landmarks,
    horizons,
    method,
    dynamic_covariates = c(),
    include_clusters = FALSE,
    validation_fold = 0,
    ...
  ) {
    # Check that method is a function with arguments formula, data, ...
    if (is(method)[1] == "character" && method == "coxph") {
      method <- predict
    }
    if (!(is(method)[1] == "function")) {
      stop(
        "Argument method",
        " must be a function",
        "\n"
      )
    }
    if (length(landmarks) != length(horizons)) {
      stop("@landmarks and @horizons must be of the same length.")
    } else if (length(landmarks) == 1) {
      # Base case for recursion
      # Check that relevant risk set is available
      if (!(landmarks %in% x@landmarks)) {
        stop(
          "Risk set for landmark time ",
          landmarks,
          " has not been computed\n"
        )
      }
      model_name <- paste0(landmarks, "-", horizons)
      # Check that relevant model fit is available
      if (!(model_name %in% names(x@survival_fits))) {
        stop(
          "Survival model has not been fitted for horizons",
          horizons,
          " at landmark time",
          landmarks,
          "\n"
        )
      }

      x@survival_predictions[[model_name]] <- method(
        x@survival_fits[[model_name]],
        ...
      )
      # Out-of-sample predictions
      if (validation_fold > 0) {
        x@survival_datasets_test[[model_name]] <-
          .create_survival_dataframe(
            x,
            landmarks,
            horizons,
            dynamic_covariates,
            include_clusters,
            validation_fold,
            train = FALSE
          )
        x@survival_predictions_test[[model_name]] <- method(
          x@survival_fits[[model_name]],
          x@survival_datasets_test[[model_name]],
          ...
        )
      }
    } else {
      # Recursion
      x <- predict_survival(
        x,
        landmarks[1],
        horizons[1],
        method,
        dynamic_covariates,
        include_clusters,
        validation_fold,
        ...
      )
      x <- predict_survival(
        x,
        landmarks[-1],
        horizons[-1],
        method,
        dynamic_covariates,
        include_clusters,
        validation_fold,
        ...
      )
    }
    x
  }
)
