#' Fits the specified longitudinal model for time-varying covariates up to
#' the landmark times
#'
#' @details
#'   ## Parallel processing
#'   As the longitudinal model for each landmark time is independent of
#'   the longitudinal models for other landmark times, parallel processing can
#'   be used to vastly speed up computation. However, due to issues with
#'   parallel processing in R, currently only Unix-like operating systems
#'   are supported by \code{landmaRk}.
#'
#' @param x An object of class \code{\link{Landmarking}}.
#' @param landmarks A vector of Landmark times.
#' @param method Either \code{"lcmm"} or \code{"lme4"} or a function for fitting
#'   a longitudinal data model, where the first argument is a formula, and also
#'   has a \code{data} argument.
#' @param formula A formula to be used in longitudinal sub-model fitting.
#' @param dynamic_covariates Vector of time-varying covariates to be modelled
#'   as the outcome of a longitudinal model.
#' @param cores Number of cores/threads to be used for parallel computation on
#'   Linux and MacOS. Defaults to either \code{options("Ncpus")} if set, or 1
#'   (single threaded) otherwise. Only single-threaded computation is currently
#'   supported on Windows.
#' @param ... Additional arguments passed to the longitudinal model fitting
#'   function (e.g. number of classes/clusters for lcmm).
#' @returns An object of class \code{\link{Landmarking}}.
#' @export
#' @seealso [lcmm::hlme()] and [lme4::lmer()] for additional arguments.
#' @examples
setGeneric(
  "fit_longitudinal",
  function(
    x,
    landmarks,
    method,
    formula,
    dynamic_covariates,
    cores = getOption("Ncpus", 1L),
    ...
  ) {
    standardGeneric("fit_longitudinal")
  }
)

#' Fits the specified longitudinal model for time-varying covariates up to
#' the landmark times
#'
#' @details
#'   ## Parallel processing
#'   As the longitudinal model for each landmark time is independent of
#'   the longitudinal models for other landmark times, parallel processing can
#'   be used to vastly speed up computation. However, due to issues with
#'   parallel processing in R, currently only Unix-like operating systems
#'   are supported by \code{landmaRk}.
#'
#' @inheritParams fit_longitudinal
#' @returns An object of class \code{\link{Landmarking}}.
#' @seealso [lcmm::hlme()] and [lme4::lmer()] for additional arguments.
#' @export
#'
#' @examples
setMethod(
  "fit_longitudinal",
  "Landmarking",
  function(
    x,
    landmarks,
    method,
    formula,
    dynamic_covariates,
    cores = getOption("Ncpus", 1L),
    ...
  ) {
    landmark <- NULL # Global var

    method <- check_method_long_fit(method)

    if (Sys.info()["sysname"] == "Windows") {
      `%doparallel%` <- foreach::`%do%`
    } else {
      cl <- init_cl(cores)
      `%doparallel%` <- foreach::`%dopar%`
      on.exit(parallel::stopCluster(cl), add = TRUE)
    }

    # If more than 25% of the individuals have 0 or 1 observations,
    # raise a warning
    for (landmark in landmarks) {
      for (dynamic_covariate in dynamic_covariates) {
        at_risk_individuals <- x@risk_sets[[as.character(landmark)]]
        # Construct dataset for the longitudinal analysis (static measurements +
        # time-varying covariate and its recording time)
        dataframe <- construct_data(
          x,
          dynamic_covariate,
          at_risk_individuals,
          landmark
        )
        prop_individuals_few_obs <- sum(table(dataframe$id) <= 1) /
          length(at_risk_individuals)
        if (prop_individuals_few_obs >= 0.25) {
          warning(
            round(prop_individuals_few_obs * 100, 2),
            "% of the individuals have 0",
            " or 1 observations at landmark time ",
            landmark,
            " for ",
            "longitudinal covariate ",
            dynamic_covariate
          )
        }
      }
    }

    x@longitudinal_fits <- foreach::foreach(landmark = landmarks) %doparallel%
      {
        check_riskset(x, landmark)
        # Create list for storing model fits for longitudinal analysis
        model_fits <- list()

        # Risk set for the landmark time
        at_risk_individuals <- x@risk_sets[[as.character(landmark)]]
        # Loop that iterates over all time-varying covariates to fit a longitudinal
        # model for the underlying trajectories
        for (dynamic_covariate in dynamic_covariates) {
          check_dynamic_covariate(x, dynamic_covariate)
          # Construct dataset for the longitudinal analysis (static measurements +
          # time-varying covariate and its recording time)
          dataframe <- construct_data(
            x,
            dynamic_covariate,
            at_risk_individuals,
            landmark
          )

          # Fit longitudinal model according to chosen method
          model_fits[[
            dynamic_covariate
          ]] <- method(
            formula,
            data = dataframe,
            ...
          )
        }
        model_fits
      }
    names(x@longitudinal_fits) <- landmarks
    x
  }
)


#' Make predictions for time-varying covariates at specified landmark times
#'
#' @param x An object of class \code{\link{Landmarking}}.
#' @param landmarks A numeric vector of landmark times.
#' @param method Longitudinal data analysis method used to make predictions
#' @param dynamic_covariates Vector of time-varying covariates to be modelled
#'   as the outcome of a longitudinal model.
#' @param ... Additional arguments passed to the prediction function (e.g.
#'   number of classes/clusters for lcmm).
#'
#' @returns An object of class \code{\link{Landmarking}}.
#' @export
#'
#' @examples
setGeneric(
  "predict_longitudinal",
  function(x, landmarks, method, dynamic_covariates, ...) {
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
  "predict_longitudinal",
  "Landmarking",
  function(x, landmarks, method, dynamic_covariates, ...) {
    value <- NULL # Global var

    method <- check_method_long_predict(method)

    # Base case for recursion
    if (length(landmarks) == 1) {
      # Check that relevant risk set is available
      check_riskset(x, landmarks)
      if (inherits(method, "character") && method == "locf") {
        # If model method is LOCF
        # Relevant risk set
        risk_set <- x@risk_sets[[as.character(landmarks)]]
        for (dynamic_covariate in dynamic_covariates) {
          predictions <- as.data.frame(risk_set)
          colnames(predictions) <- x@ids
          predictions <- predictions |>
            left_join(
              x@data_dynamic[[dynamic_covariate]] |>
                filter(get(x@times) <= landmarks) |>
                slice_max(get(x@times), by = id),
              by = stats::setNames(x@ids, x@ids)
            )
          predictions <- predictions |> pull(x@measurements, name = x@ids)
          # Impute NAs
          if (any(is.na(predictions))) {
            if (is.numeric(predictions)) {
              # Replace NAs with mean if covariate is continuous
              predictions[is.na(predictions)] <- mean(predictions, na.rm = TRUE)
            } else {
              # Replace NAs with mode if covariate is discrete
              predictions[is.na(predictions)] <- names(sort(
                -table(predictions)
              ))[1]
            }
          }
          x@longitudinal_predictions[[as.character(landmarks)]][[
            dynamic_covariate
          ]] <- predictions
        }
      } else {
        # Check that relevant model fit is available (if model is not LOCF)
        check_long_fit(x, landmarks)

        # Relevant risk set
        risk_set <- x@risk_sets[[as.character(landmarks)]]
        # Create list for storing model predictions, for longitudinal analysis
        x@longitudinal_predictions[[as.character(landmarks)]] <- list()
        # Loop that iterates over all time-varying covariates, to fit a longitudinal
        # model for the underlying trajectories
        for (dynamic_covariate in names(x@data_dynamic)) {
          # Check that relevant model fit is available
          if (
            !(dynamic_covariate %in%
              names(x@longitudinal_fits[[as.character(landmarks)]]))
          ) {
            stop(
              "Longitudinal model has not been fit for dynamic covariate ",
              dynamic_covariate,
              " at landmark time",
              landmarks,
              ". Fit a longitudinal model using fit_longitudinal,",
              "or use Last Observation Carried Forward (LOCF).",
              "\n"
            )
          } else {
            # Fit longitudinal model according to chosen method
            newdata <- data.frame(risk_set, landmarks)
            colnames(newdata) <- c(x@ids, x@times)
            newdata <- newdata |>
              left_join(x@data_static, by = stats::setNames(x@ids, x@ids))

            x@longitudinal_predictions[[as.character(landmarks)]][[
              dynamic_covariate
            ]] <- method(
              x@longitudinal_fits[[as.character(landmarks)]][[
                dynamic_covariate
              ]],
              newdata = newdata,
              ...
            )
          }

          x@longitudinal_predictions[[as.character(landmarks)]][[
            dynamic_covariate
          ]] <- method(
            x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
            newdata = newdata,
            ...
          )
          if (
            length(x@longitudinal_predictions[[as.character(landmarks)]][[
              dynamic_covariate
            ]]) !=
            nrow(newdata)
          ) {
            stop(paste(
              "Number of predictions for dynamic_covariate",
              dynamic_covariate,
              "at landmark time",
              landmarks,
              "differs from number of observations in the risk set."
            ))
          }
        }
      }
    } else {
      # Recursion
      x <- predict_longitudinal(
        x,
        landmarks[1],
        method,
        dynamic_covariates,
        ...
      )
      x <- predict_longitudinal(
        x,
        landmarks[-1],
        method,
        dynamic_covariates,
        ...
      )
    }
    x
  }
)
