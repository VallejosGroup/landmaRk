# Check that method is a function with arguments formula, data, ...
.check_method_long_fit <- function(method) {
  if (is(method)[1] == "character" && method == "lcmm") {
    method <- .fit_lcmm
  } else if (is(method)[1] == "character" && method == "lme4") {
    method <- lme4::lmer
  }

  if (is(method)[1] != "function") {
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
  method
}

.check_method_long_predict <- function(method) {
  # Check that method is a function with arguments formula, data, ...
  if (is(method)[1] == "character" && method == "lcmm") {
    method <- .predict_lcmm
  } else if (is(method)[1] == "character" && method == "lme4") {
    method <- .predict_lme4
  } else if (is(method)[1] == "character" && method == "locf") {
    method <- .locf_summary
  } else if (!(is(method)[1] == "function")) {
    stop(
      "Argument method must be one of 'lme4', 'lcmm',",
      " 'locf' or a function",
      "\n"
    )
  }
  method
}

# Check whether `method` is a "summary measure": a function that computes a
# longitudinal summary directly from the raw data (e.g. LOCF, or a
# user-supplied last value/mean/etc. function), without requiring a model
# to have been previously fit with fit_longitudinal(). Such functions are
# identified by having `data`, `id`, `time`, `value` and `landmark` among
# their formal arguments, as opposed to fit-based prediction functions
# (e.g. lcmm, lme4) which take a fitted model object as their first
# argument.
.is_summary_method <- function(method) {
  is.function(method) &&
    all(
      c("data", "id", "time", "value", "landmark") %in% names(formals(method))
    )
}

# Check that a summary-measure method is a function with the arguments
# required by .compute_summary_predictions(): data, id, time, value,
# landmark, ...
.check_method_long_summary <- function(method) {
  if (!is.function(method)) {
    stop(
      "Argument ",
      method,
      " must be a function",
      "\n"
    )
  }
  required_args <- c("data", "id", "time", "value", "landmark")
  missing_args <- setdiff(required_args, names(formals(method)))
  if (length(missing_args) > 0) {
    stop(
      "A summary-measure method must be a function with arguments ",
      paste(required_args, collapse = ", "),
      ". Missing: ",
      paste(missing_args, collapse = ", "),
      "\n"
    )
  }
  invisible(method)
}


.check_riskset <- function(x, landmark) {
  # Check that relevant risk set is available
  if (!(landmark %in% x@landmarks)) {
    stop(
      "Risk set for landmark time ",
      landmark,
      " has not been computed",
      "\n"
    )
  }
}

.check_dynamic_covariate <- function(x, dynamic_covariate) {
  if (!(dynamic_covariate %in% names(x@data_dynamic))) {
    stop(
      "Data frame has not been provided for dynamic covariate",
      dynamic_covariate
    )
  }
}

# Check that longitudinal model is available for prediction
.check_long_fit <- function(x, landmarks) {
  # Check that relevant model fit is available
  if (!(as.character(landmarks) %in% names(x@longitudinal_fits))) {
    stop(
      "Longitudinal model has not been fit for landmark time",
      landmarks,
      "\n"
    )
  }
}

# Construct data frame for longitudinal model fitting
.construct_data <- function(
  x,
  dynamic_covariate,
  at_risk_individuals,
  landmark
) {
  if (x@censor_at_landmark) {
    at_risk_individuals <- data.frame(at_risk_individuals)
    colnames(at_risk_individuals) <- x@ids
    if (inherits(x@data_dynamic[[dynamic_covariate]], "tbl_df")) {
      at_risk_individuals <- dplyr::as_tibble(at_risk_individuals)
    }
  } else {
    at_risk_individuals <- data.frame(x@data_static[, x@ids])
    colnames(at_risk_individuals) <- x@ids
  }

  result_df <- at_risk_individuals |>
    # Subset with individuals who are at risk only
    left_join(
      x@data_dynamic[[dynamic_covariate]],
      by = stats::setNames(x@ids, x@ids)
    ) |>
    # Join with static covariates
    dplyr::left_join(x@data_static, by = x@ids)

  if (x@censor_at_landmark) {
    # Subset with observations prior to landmark time
    result_df <- result_df |>
      dplyr::filter(get(x@times) <= landmark)
  }
  result_df
}

# Initialize a cluster for parallel processing based on the operating system
.init_cl <- function(cores) {
  # Use FORK on Unix-like systems
  cl <- parallel::makeCluster(cores, type = "FORK")
  doParallel::registerDoParallel(cl)
  cl
}

# Default summary measure: Last Observation Carried Forward. Picks, for
# each individual, the most recent measurement recorded at or before the
# landmark time. See .compute_summary_predictions() for the calling
# convention shared by all summary-measure methods.
.locf_summary <- function(data, id, time, value, landmark, ...) {
  data |>
    dplyr::slice_max(get(time), by = dplyr::all_of(id)) |>
    dplyr::pull(value, name = id)
}

# Compute summary-measure predictions for a given fold condition, using an
# arbitrary user-supplied summary function (e.g. LOCF, or a custom
# mean/last-value/etc. function). Unlike fit-based methods (lcmm, lme4),
# summary methods operate directly on the raw longitudinal data and do not
# require a model to have been previously fit with fit_longitudinal() (see
# .is_summary_method() / .check_method_long_summary()).
.compute_summary_predictions <- function(
  x,
  risk_set,
  dynamic_covariate,
  landmarks,
  cv_folds_subset,
  method,
  ...
) {
  .check_method_long_summary(method)

  at_risk_ids <- as.data.frame(risk_set)
  colnames(at_risk_ids) <- x@ids
  at_risk_ids <- at_risk_ids |>
    dplyr::inner_join(
      cv_folds_subset |> dplyr::select(x@ids),
      by = x@ids
    )

  # Longitudinal data for the at-risk individuals, observed at or before the
  # landmark time
  data <- x@data_dynamic[[dynamic_covariate]] |>
    dplyr::filter(get(x@times) <= landmarks) |>
    dplyr::inner_join(at_risk_ids, by = x@ids)

  predictions <- method(
    data = data,
    id = x@ids,
    time = x@times,
    value = x@measurements,
    landmark = landmarks,
    ...
  )

  if (is.data.frame(predictions)) {
    predictions <- predictions |> dplyr::pull(x@measurements, name = x@ids)
  }

  # Re-index to the full set of at-risk individuals: those for which
  # `method` did not return a prediction (e.g. no observations before the
  # landmark time) become NA, and are imputed below
  ids <- at_risk_ids[[x@ids]]
  predictions <- predictions[as.character(ids)]
  names(predictions) <- ids

  # Impute NAs
  if (any(is.na(predictions))) {
    n_imputed <- sum(is.na(predictions))
    message(
      n_imputed,
      " individual(s) have no observations before landmark time ",
      landmarks,
      " for dynamic covariate ",
      dynamic_covariate,
      ". Imputing with population ",
      if (is.numeric(predictions)) "mean." else "mode."
    )
    if (is.numeric(predictions)) {
      # Replace NAs with mean if covariate is continuous
      predictions[is.na(predictions)] <- mean(predictions, na.rm = TRUE)
    } else {
      # Replace NAs with mode if covariate is discrete
      predictions[is.na(predictions)] <- names(sort(-table(predictions)))[1]
    }
  }
  predictions
}

# Check if system supports parallel processing
.supports_parallel <- function() {
  Sys.info()[["sysname"]] != "Windows"
}

.fit_longitudinal_model <- function(
  x,
  landmark,
  method,
  formula,
  dynamic_covariates,
  validation_fold = 0,
  ...
) {
  fold <- NULL
  .check_riskset(x, landmark)
  # Create list for storing model fits for longitudinal analysis
  model_fits <- list()

  # Risk set for the landmark time
  at_risk_individuals <- x@risk_sets[[as.character(landmark)]]
  # Loop that iterates over all time-varying covariates to fit a
  # longitudinal model for the underlying trajectories
  for (dynamic_covariate in dynamic_covariates) {
    .check_dynamic_covariate(x, dynamic_covariate)
    # Construct dataset for the longitudinal analysis
    # (static measurements+time-varying covariate and its recording time)
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
