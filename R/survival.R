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
#' @param x An object of class \code{\link{Landmarking}}.
#' @param formula A formula to be used in survival sub-model fitting.
#' @param landmarks Numeric vector of landmark times.
#' @param horizons Vector of prediction horizons up to when the survival submodel
#'   is fitted.
#' @param method Method for survival analysis, either "survfit" or "coxph".
#' @param dynamic_covariates Vector of time-varying covariates to be used
#' in the survival model.
#'
#' @returns An object of class Landmarking.
#' @export
#'
#' @examples
setGeneric(
  "fit_survival",
  function(x, formula, landmarks, horizons, method, dynamic_covariates = c()) {
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
#' @returns An object of class Landmarking.
#' @export
#'
#' @examples
setMethod(
  "fit_survival",
  "Landmarking",
  function(x, formula, landmarks, horizons, method, dynamic_covariates = c()) {
    # Check that method is a function with arguments formula, data, ...
    method <- .check_method_survival_predict(method)

    if (length(landmarks) != length(horizons)) {
      stop("@landmarks and @horizons must be of the same length.")
    } else if (length(landmarks) == 1) {
      # Base case for recursion
      .check_riskset_survival(x, landmarks)
      # Recover risk sets (ids of individuals who are at risk at landmark time)
      at_risk_individuals <- x@risk_sets[[as.character(landmarks)]]

        # Construct dataset for survival analysis (censor events past horizon time)
        x@survival_datasets[[paste0(landmarks, "-", horizons)]] <- x@data_static[
          which(x@data_static[, x@event_time] >= landmarks),
        ] |>
          mutate(
            event_status = ifelse(
              get(x@event_time) > horizons,
              0,
              get(x@event_indicator)
            ),
            event_time = ifelse(
              get(x@event_time) > horizons,
              horizons - landmarks,
              get(x@event_time) - landmarks
            )
          )

        # Check that longitudinal predictions are available at landmark time
        if (length(dynamic_covariates) > 0) {
          .check_predictions_available_survival(
            x,
            landmarks,
            dynamic_covariates
          )

          x@survival_datasets[[paste0(landmarks, "-", horizons)]] <- cbind(
            x@survival_datasets[[paste0(landmarks, "-", horizons)]],
            do.call(
              bind_cols,
              x@longitudinal_predictions[[as.character(landmarks)]]
            )
          )
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
        dynamic_covariates
      )
      x <- fit_survival(
        x,
        formula,
        landmarks[-1],
        horizons[-1],
        method,
        dynamic_covariates
      )
    }
    x
  }
)

#' Make predictions for time-to-event outcomes at specified horizon times
#'
#' @param x An object of class \code{\link{Landmarking}}.
#' @param landmarks A numeric vector of landmark times.
#' @param horizons Vector of prediction horizons up to when the survival submodel
#'   is fitted.
#' @param method R function that is used to make predictions
#' @param ... Additional arguments passed to the prediction function (e.g.
#'   number of classes/clusters for lcmm).
#'
#' @returns An object of class \code{\link{Landmarking}}.
#' @export
#'
#' @examples
setGeneric(
  "predict_survival",
  function(x, landmarks, horizons, method, ...) {
    standardGeneric("predict_survival")
  }
)

#' Make predictions for time-to-event outcomes at specified horizon times
#'
#' @inheritParams predict_survival
#'
#' @returns An object of class \code{\link{Landmarking}}.
#' @export
#'
#' @examples
setMethod(
  "predict_survival",
  "Landmarking",
  function(x, landmarks, horizons, method, ...) {
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
    } else {
      # Recursion
      x <- predict_survival(x, landmarks[1], horizons[1], method, ...)
      x <- predict_survival(x, landmarks[-1], horizons[-1], method, ...)
    }
    x
  }
)
