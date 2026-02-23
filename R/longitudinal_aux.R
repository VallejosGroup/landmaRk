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
    method <- "locf"
  } else if (!(is(method)[1] == "function")) {
    stop(
      "Argument method must be one of 'lme4', 'lcmm',",
      " 'locf' or a function",
      "\n"
    )
  }
  method
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
  censor_at_landmark,
  landmark
) {
  if (censor_at_landmark) {
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

  if (censor_at_landmark) {
    # Subset with observations prior to landmark time
    result_df <- result_df |>
      dplyr::filter(get(x@times) <= landmark)
  }
  return(result_df)
}

# Initialize a cluster for parallel processing based on the operating system
.init_cl <- function(cores) {
  # Use FORK on Unix-like systems
  cl <- parallel::makeCluster(cores, type = "FORK")
  doParallel::registerDoParallel(cl)
  cl
}

# Compute LOCF predictions for a given fold condition
.compute_locf_predictions <- function(
  x,
  risk_set,
  dynamic_covariate,
  landmarks,
  cv_folds_subset
) {
  predictions <- as.data.frame(risk_set)
  colnames(predictions) <- x@ids
  predictions <- predictions |>
    dplyr::inner_join(
      cv_folds_subset |> dplyr::select(x@ids),
      by = x@ids
    ) |>
    dplyr::left_join(
      x@data_dynamic[[dynamic_covariate]] |>
        dplyr::filter(get(x@times) <= landmarks) |>
        dplyr::slice_max(get(x@times), by = x@ids),
      by = stats::setNames(x@ids, x@ids)
    )
  predictions <- predictions |> dplyr::pull(x@measurements, name = x@ids)
  # Impute NAs
  if (any(is.na(predictions))) {
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
