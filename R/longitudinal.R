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
#' @param x An object of class \code{\link{LandmarkAnalysis}}.
#' @param landmarks A vector of Landmark times.
#' @param method Either \code{"lcmm"} or \code{"lme4"} or a function for fitting
#'   a longitudinal data model, where the first argument is a formula, and also
#'   has a \code{data} argument. Only needed for fit-based prediction methods
#'   used later in \code{\link{predict_longitudinal}} (e.g. \code{"lcmm"},
#'   \code{"lme4"}); summary measures such as \code{"locf"} are computed
#'   directly from the data and do not require a call to
#'   \code{fit_longitudinal}.
#' @param formula A formula to be used in longitudinal sub-model fitting.
#' @param dynamic_covariates Vector of time-varying covariates to be modelled
#'   as the outcome of a longitudinal model.
#' @param validation_fold If positive, cross-validation fold where model is
#'   fitted. If 0 (default), model fitting is performed using the complete
#'   dataset.
#' @param cores Number of cores/threads to be used for parallel computation on
#'   Linux and MacOS. Defaults to either \code{options("Ncpus")} if set, or 1
#'   (single threaded) otherwise. Only single-threaded computation is currently
#'   supported on Windows.
#' @param .warn_when_prop_few_obs Threshold proportion (0-1) for warning when
#'   individuals have 0 or 1 observations. Defaults to 0.25 (i.e., warn when
#'   25% or more individuals have few observations).
#' @param ... Additional arguments passed to the longitudinal model fitting
#'   function (e.g. number of classes/clusters for lcmm).
#' @returns An object of class \code{\link{LandmarkAnalysis}}.
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
    validation_fold = 0,
    cores = getOption("Ncpus", 1L),
    .warn_when_prop_few_obs = 0.25,
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
#' @returns An object of class \code{\link{LandmarkAnalysis}}.
#' @seealso [lcmm::hlme()] and [lme4::lmer()] for additional arguments.
#' @export
#'
#' @examples
setMethod(
  "fit_longitudinal",
  "LandmarkAnalysis",
  function(
    x,
    landmarks,
    method,
    formula,
    dynamic_covariates,
    validation_fold = 0,
    cores = getOption("Ncpus", 1L),
    .warn_when_prop_few_obs = 0.25,
    ...
  ) {
    landmark <- NULL # Global var
    fold <- NULL # Global var

    if (
      !is.numeric(.warn_when_prop_few_obs) ||
        length(.warn_when_prop_few_obs) != 1L ||
        is.na(.warn_when_prop_few_obs) ||
        .warn_when_prop_few_obs < 0 ||
        .warn_when_prop_few_obs > 1
    ) {
      stop(
        "@.warn_when_prop_few_obs must be a single numeric value between 0 and 1"
      )
    }

    method <- .check_method_long_fit(method)

    if (.supports_parallel()) {
      cl <- .init_cl(cores)
      `%doparallel%` <- foreach::`%dopar%`
      on.exit(parallel::stopCluster(cl), add = TRUE)
    } else {
      `%doparallel%` <- foreach::`%do%`
    }

    # Warn if a large proportion of individuals have 0 or 1 observations
    for (landmark in landmarks) {
      for (dynamic_covariate in dynamic_covariates) {
        at_risk_individuals <- x@risk_sets[[as.character(landmark)]]
        # Construct dataset for the longitudinal analysis (static measurements +
        # time-varying covariate and its recording time)
        dataframe <- .construct_data(
          x,
          dynamic_covariate,
          at_risk_individuals,
          landmark
        ) |>
          inner_join(
            x@cv_folds |> filter(fold != validation_fold) |> select(x@ids),
            by = x@ids
          )
        prop_individuals_few_obs <- sum(table(dataframe[, x@ids]) <= 1) /
          length(at_risk_individuals)
        if (prop_individuals_few_obs >= .warn_when_prop_few_obs) {
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

    if (!x@censor_at_landmark) {
      longitudinal_fits <- .fit_longitudinal_model(
        x,
        landmarks[1],
        method,
        formula,
        dynamic_covariates,
        validation_fold,
        ...
      )

      x@longitudinal_fits <- lapply(as.character(landmarks), function(i) {
        longitudinal_fits
      })
    } else {
      x@longitudinal_fits <-
        foreach::foreach(landmark = landmarks) %doparallel%
        {
          .fit_longitudinal_model(
            x,
            landmark,
            method,
            formula,
            dynamic_covariates,
            validation_fold,
            ...
          )
        }
    }
    names(x@longitudinal_fits) <- landmarks
    x
  }
)


#' Make predictions for time-varying covariates at specified landmark times
#'
#' @param x An object of class \code{\link{LandmarkAnalysis}}.
#' @param landmarks A numeric vector of landmark times.
#' @param method Longitudinal data analysis method used to make
#'   predictions. Either \code{"lcmm"}, \code{"lme4"}, \code{"locf"}, or a
#'   function, which can be one of two kinds:
#'   \itemize{
#'     \item A summary measure, like \code{"locf"}, computed directly from
#'       the raw longitudinal data and not requiring a model to have been
#'       previously fit with \code{\link{fit_longitudinal}}. Such a
#'       function must have the arguments \code{data}, \code{id},
#'       \code{time}, \code{value} and \code{landmark} (and optionally
#'       further arguments passed through \code{...}), and must return a
#'       named vector (or a two-column data frame) with one summary value
#'       per individual in the risk set.
#'     \item A prediction function for a model previously fit with
#'       \code{\link{fit_longitudinal}} (as is the case for \code{"lcmm"}
#'       and \code{"lme4"}), where the first argument is the fitted model
#'       object, and which also has \code{newdata} and \code{subject}
#'       arguments.
#'   }
#' @param dynamic_covariates Vector of time-varying covariates to be modelled
#'   as the outcome of a longitudinal model.
#' @param validation_fold If positive, cross-validation fold where model is
#'   fitted. If 0 (default), model fitting is performed in the complete dataset.
#' @param ... Additional arguments passed to the prediction function (e.g.
#'   number of classes/clusters for lcmm).
#'
#' @returns An object of class \code{\link{LandmarkAnalysis}}.
#' @export
#'
#' @examples
setGeneric(
  "predict_longitudinal",
  function(
    x,
    landmarks,
    method,
    dynamic_covariates,
    validation_fold = 0,
    ...
  ) {
    standardGeneric("predict_longitudinal")
  }
)

#' Make predictions for time-varying covariates at specified landmark times
#'
#' @inheritParams predict_longitudinal
#'
#' @returns An object of class \code{\link{LandmarkAnalysis}}.
#' @export
#'
#' @examples
setMethod(
  "predict_longitudinal",
  "LandmarkAnalysis",
  function(
    x,
    landmarks,
    method,
    dynamic_covariates,
    validation_fold = 0,
    ...
  ) {
    fold <- NULL # Global var

    method <- .check_method_long_predict(method)

    # Base case for recursion
    if (length(landmarks) == 1) {
      # Check that relevant risk set is available
      .check_riskset(x, landmarks)
      if (.is_summary_method(method)) {
        # If method is a summary measure (e.g. "locf", or a custom
        # function), computed directly from the raw longitudinal data, with
        # no model fit required
        # Relevant risk set
        risk_set <- x@risk_sets[[as.character(landmarks)]]
        for (dynamic_covariate in dynamic_covariates) {
          # Training fold predictions
          train_folds <- x@cv_folds |> filter(fold != validation_fold)
          x@longitudinal_predictions[[as.character(landmarks)]][[
            dynamic_covariate
          ]] <- .compute_summary_predictions(
            x,
            risk_set,
            dynamic_covariate,
            landmarks,
            train_folds,
            method,
            ...
          )

          # Test fold predictions (if cross-validation is enabled)
          if (validation_fold > 0) {
            test_folds <- x@cv_folds |> filter(fold == validation_fold)
            x@longitudinal_predictions_test[[as.character(landmarks)]][[
              dynamic_covariate
            ]] <- .compute_summary_predictions(
              x,
              risk_set,
              dynamic_covariate,
              landmarks,
              test_folds,
              method,
              ...
            )
          }
        }
      } else {
        # Check that relevant model fit is available (only required for
        # fit-based methods, e.g. "lcmm"/"lme4" or a custom predict-from-fit
        # function; summary measures are handled above and need no prior
        # fit_longitudinal() call)
        .check_long_fit(x, landmarks)

        # Relevant risk set
        risk_set <- x@risk_sets[[as.character(landmarks)]]
        # Create list for storing model predictions, for longitudinal analysis
        x@longitudinal_predictions[[as.character(landmarks)]] <- list()
        if (validation_fold > 0) {
          x@longitudinal_predictions_test[[as.character(landmarks)]] <- list()
        }
        # Loop that iterates over all time-varying covariates, to fit a
        # longitudinal model for the underlying trajectories
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

            newdata_train <- newdata |>
              inner_join(
                x@cv_folds |> filter(fold != validation_fold) |> select(x@ids),
                by = x@ids
              )
            x@longitudinal_predictions[[as.character(landmarks)]][[
              dynamic_covariate
            ]] <- method(
              x@longitudinal_fits[[as.character(landmarks)]][[
                dynamic_covariate
              ]],
              newdata = newdata_train,
              subject = x@ids,
              ...
            )

            if (validation_fold > 0) {
              newdata <- newdata |>
                inner_join(
                  x@cv_folds |>
                    filter(fold == validation_fold) |>
                    select(x@ids),
                  by = x@ids
                )
              # Check optional arguments
              mc <- match.call(expand.dots = FALSE)$...
              # Vector nms contains
              nms <- names(mc)
              if (
                !is.null(nms) &&
                  (("include_clusters" %in%
                    names(list(...)) &&
                    list(...)$include_clusters) ||
                    ("avg" %in% names(list(...)) && list(...)$avg))
              ) {
                newdata <- newdata |>
                  left_join(
                    x@data_dynamic[[dynamic_covariate]] |>
                      filter(get(x@times) <= landmarks) |>
                      select(-any_of(x@times)) |>
                      select(-any_of(x@measurements)) |>
                      unique(),
                    by = x@ids
                  )
              }
              x@longitudinal_predictions_test[[as.character(landmarks)]][[
                dynamic_covariate
              ]] <- method(
                x@longitudinal_fits[[as.character(landmarks)]][[
                  dynamic_covariate
                ]],
                newdata = newdata |>
                  dplyr::left_join(
                    x@data_dynamic[[dynamic_covariate]] |>
                      dplyr::filter(get(x@times) <= landmarks) |>
                      dplyr::slice_max(get(x@times), by = x@ids) |>
                      select(-!!sym(x@times)),
                    by = stats::setNames(x@ids, x@ids)
                  ),
                subject = x@ids,
                test = TRUE,
                newdata_long = newdata |>
                  select(-any_of(x@times)) |>
                  select(-any_of(x@measurements)) |>
                  inner_join(
                    x@data_dynamic[[dynamic_covariate]] |>
                      filter(get(x@times) <= landmarks),
                    by = x@ids
                  ),
                ...
              )
            }
          }

          predictions <- x@longitudinal_predictions[[as.character(landmarks)]][[
            dynamic_covariate
          ]]
          # Number of predictions (length if stored in vector or number of rows
          # if stored in matrix)
          npred <- ifelse(
            is.null(dim(predictions)),
            length(predictions),
            nrow(predictions)
          )
          if (npred != nrow(newdata_train)) {
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
        validation_fold,
        ...
      )
      x <- predict_longitudinal(
        x,
        landmarks[-1],
        method,
        dynamic_covariates,
        validation_fold,
        ...
      )
    }
    x
  }
)
