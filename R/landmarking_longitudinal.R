#' Fits the specified longitudinal model for the latent processes underlying the
#' relevant time-varying covariates, up until the landmarking times
#'
#' @param x An object of class \code{\link{Landmarking}}.
#' @param landmarks A vector of Landmark times.
#' @param method A function for fitting a longitudinal data model, whose first
#' argument is a formula, and has at least an argument "data".
#' @param formula A formula to be used in longitudinal submodel fitting.
#' @param ... Additional arguments passed to the longitudinal model fitting
#'   function (e.g. number of classes/clusters for lcmm).
#' @returns An object of class \code{\link{Landmarking}}.
#' @export
#'
#' @examples
setGeneric(
  "fit_longitudinal",
  function(x,
           landmarks,
           method,
           formula,
           ...) {
    standardGeneric("fit_longitudinal")
  }
)

#' Fits the specified longitudinal model for the latent processes underlying the
#' relevant time-varying covariates, up until the landmarking times
#'
#' @inheritParams fit_longitudinal
#' @returns An object of class \code{\link{Landmarking}}.
#' @export
#'
#' @examples
setMethod(
  "fit_longitudinal",
  "Landmarking",
  function(x, landmarks, method, formula, ...) {
    # Check that method is a function with arguments formula, data, ...
    if (is(method)[1] == "character" && method == "lcmm") {
      method <- fit_lcmm_
    } else if (is(method)[1] == "character" && method == "lme4") {
      method <- lme4::lmer
    }
    if (!(is(method)[1] == "function")) {
      stop(
        "Argument ",
        method,
        " must be a function",
        "\n"
      )
    }
    if (!("data" %in% names(as.list(args(method))))) {
      stop(
        "Argument ",
        method,
        " must be a function, and data must be an argument to that function",
        "\n"
      )
    }
    # Base case for recursion
    if (length(landmarks) == 1) {
      # Check that relevant risk set is available
      if (!(landmarks %in% x@landmarks)) {
        stop(
          "Risk set for landmark time ",
          landmarks,
          " has not been computed",
          "\n"
        )
      }
      # Create list for storing model fits for longitudinal analysis
      x@longitudinal_fits[[as.character(landmarks)]] <- list()

      # Risk set for the landmark time
      at_risk_individuals <- x@risk_sets[[as.character(landmarks)]]
      # Loop that iterates over all time-varying covariates to fit a longitudinal
      # model for the underlying trajectories
      for (predictor in names(x@data_dynamic)) {
        # Construct dataset for the longitudinal analysis (static measurements +
        # time-varying covariate and its recording time)
        dataframe <- x@data_dynamic[[predictor]] |>
          # Subset with individuals who are at risk only
          filter(get(x@ids) %in% at_risk_individuals) |>
          # Subset with observations prior to landmark time
          filter(get(x@times) <= landmarks) |>
          # Join with static covariates
          left_join(x@data_static, by = x@ids)
        # Fit longitudinal model according to chosen method
        x@longitudinal_fits[[as.character(landmarks)]][[predictor]] <- method(
          formula,
          data = dataframe, ...
        )
      }
    } else {
      # Recursion
      x <- fit_longitudinal(x, landmarks[1], method, formula, ...)
      x <- fit_longitudinal(x, landmarks[-1], method, formula, ...)
    }
    x
  }
)

#' Make predictions for time-varying covariates at specified landmark times
#'
#' @param x An object of class \code{\link{Landmarking}}.
#' @param landmarks A numeric vector of landmark times.
#' @param method Longitudinal data analysis method used to make predictions
#' @param ... Additional arguments passed to the prediction function (e.g.
#'   number of classes/clusters for lcmm).
#'
#' @returns An object of class \code{\link{Landmarking}}.
#' @export
#'
#' @examples
setGeneric(
  "predict_longitudinal",
  function(x, landmarks, method, ...) {
    standardGeneric("predict_longitudinal")
  }
)

#' Make predictions for time-varying covariates at specified landmark times
#'
#' @inheritParams predict_longitudinal
#'
#' @returns An object of class \code{\link{Landmarking}}.
#' @export
#'
#' @examples
setMethod(
  "predict_longitudinal", "Landmarking",
  function(x, landmarks, method, ...) {
    # Check that method is a function with arguments formula, data, ...
    if (is(method)[1] == "character" && method == "lcmm") {
      method <- predict_lcmm_
    } else if (is(method)[1] == "character" && method == "lme4") {
      method <- predict
    }
    if (!(is(method)[1] == "function")) {
      stop(
        "Argument method",
        " must be a function",
        "\n"
      )
    }
    # Base case for recursion
    if (length(landmarks) == 1) {
      # Check that relevant risk set is available
      if (!(landmarks %in% x@landmarks)) {
        stop("Risk set for landmark time ", landmarks, " has not been computed\n")
      }
      # Check that relevant model fit is available
      if (!(as.character(landmarks) %in% names(x@longitudinal_fits))) {
        stop(
          "Longitudinal model has not been fit for landmark time",
          landmarks,
          "\n"
        )
      }
      # Relevant risk set
      risk_set <- x@risk_sets[[as.character(landmarks)]]
      # Create list for storing model predictions, for longitudinal analysis
      x@longitudinal_predictions[[as.character(landmarks)]] <- list()
      # Loop that iterates over all time-varying covariates, to fit a longitudinal
      # model for the underlying trajectories
      for (predictor in names(x@data_dynamic)) {
        # Check that relevant model fit is available
        if (!(predictor %in% names(x@longitudinal_fits[[as.character(landmarks)]]))) {
          stop(
            "Longitudinal model has not been fit for dynamic covariate ",
            predictor,
            " at landmark time",
            landmarks,
            "\n"
          )
        }
        # Fit longitudinal model according to chosen method
        newdata <- x@data_static |>
          filter(get(x@ids) %in% risk_set)
        newdata[, x@times] <- landmarks

        x@longitudinal_predictions[[as.character(landmarks)]][[predictor]] <- method(
          x@longitudinal_fits[[as.character(landmarks)]][[predictor]],
          newdata = newdata,
          ...
        )
        if (length(x@longitudinal_predictions[[as.character(landmarks)]][[predictor]]) != nrow(newdata)) {
          stop(paste("Number of predictions for predictor",
          predictor,
          "at landmark time",
          landmarks,
          "differs from number of observations in the risk set."))
        }
      }
    } else {
      # Recursion
      x <- predict_longitudinal(x, landmarks[1], method, ...)
      x <- predict_longitudinal(x, landmarks[-1], method, ...)
    }
    x
  }
)
