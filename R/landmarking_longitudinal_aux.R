# Check that method is a function with arguments formula, data, ...
check_method_long_fit <- function(method) {
  if (is(method)[1] == "character" && method == "lcmm") {
    method <- fit_lcmm_
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

check_method_long_predict <- function(method) {
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
  method
}


check_riskset <- function(x, landmark) {
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

check_dynamic_covariate <- function(x, dynamic_covariate) {
  if (!(dynamic_covariate %in% names(x@data_dynamic))) {
    stop(
      "Data frame has not been provided for dynamic covariate",
      dynamic_covariate
    )
  }
}

# Check that longitudinal model is available for prediction
check_long_fit <- function(x, landmarks) {
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
construct_data <- function(
  x,
  dynamic_covariate,
  at_risk_individuals,
  landmark
) {
  x@data_dynamic[[dynamic_covariate]] |>
    # Subset with individuals who are at risk only
    dplyr::filter(get(x@ids) %in% at_risk_individuals) |>
    # Subset with observations prior to landmark time
    dplyr::filter(get(x@times) <= landmark) |>
    # Join with static covariates
    dplyr::left_join(x@data_static, by = x@ids)
}


# Initialize a cluster for parallel processing based on the operating system
init_cl <- function(cores) {
  if (Sys.info()["sysname"] == "Windows") {
    # Use PSOCK on Windows
    cl <- future::makeClusterPSOCK(cores)
    doSNOW::registerDoSNOW(cl)
  } else {
    # Use FORK on Unix-like systems
    cl <- parallel::makeCluster(cores, type = "FORK")
    doParallel::registerDoParallel(cl)
  }
  cl
}
