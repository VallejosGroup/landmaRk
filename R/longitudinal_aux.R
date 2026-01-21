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
    method <- predict
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
  landmark
) {
  at_risk_individuals <- data.frame(at_risk_individuals)
  colnames(at_risk_individuals) <- x@ids
  if (inherits(x@data_dynamic[[dynamic_covariate]], "tbl_df")) {
    at_risk_individuals <- dplyr::as_tibble(at_risk_individuals)
  }

  at_risk_individuals |>
    # Subset with individuals who are at risk only
    left_join(
      x@data_dynamic[[dynamic_covariate]],
      by = stats::setNames(x@ids, x@ids)
    ) |>
    # Subset with observations prior to landmark time
    dplyr::filter(get(x@times) <= landmark) |>
    # Join with static covariates
    dplyr::left_join(x@data_static, by = x@ids)
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


# Check prediction vector length
.pred_check <- function(predRE, subject, in_train_set) {
  if (length(unique(predRE[, subject])) != length(in_train_set)) {
    stop(sprintf(
      paste0(
        "lcmm::predictRE produced %d predictions but expected %d ",
        "predictions.\n",
        "Probable reason: static covariates contain missing data.\n"
      ),
      length(unique(predRE[, subject])),
      length(in_train_set)
    ))
  }
  NULL
}


.pred_train <- function(x, hlme, subject, newdata, test) {
  predictions_step1 <- NULL

  # Step 1a. We estimate the random effects for individuals in the training set
  x$call[[1]] <- expr(hlme)
  # Step 1b. Find ids of individuals in the training set
  in_train_set <- unique(x$data[, subject])
  if (!test) {
    # In sample predictions.
    predRE <- lcmm::predictRE(x, x$data, subject = subject, classpredRE = TRUE)
    .pred_check(predRE, subject, in_train_set)

    # Step 1c. Find class-specific predictions for individuals in the
    # training set.
    predictions_step1 <- t(sapply(
      in_train_set,
      function(individual) {
        lcmm::predictY(
          x,
          newdata = newdata |> filter(get(subject) == individual),
          predRE = predRE |> filter(get(subject) == individual)
        )$pred
      }
    ))

    predictions_step1 <- as.data.frame(predictions_step1)
    predictions_step1[, subject] <- in_train_set
    predictions_step1 <- predictions_step1 |> relocate(subject)
    colnames(predictions_step1) <- c(
      subject,
      paste0("Ypred_class", 1:(ncol(predictions_step1) - 1))
    )
  }
  list(predictions_step1, in_train_set)
}


.pred_test <- function(x, subject, newdata, in_train_set) {
  predictions_step2 <- NULL
  # Step 2a. Find ids of individuals outwith the training set
  not_in_train_set <- setdiff(
    unique(newdata[, subject]),
    in_train_set
  )

  # Step 2b. Find class-specific predictions for individuals outwith the
  # training set.
  if (length(not_in_train_set) > 0) {
    predictions_step2 <- lcmm::predictY(
      x,
      newdata = newdata |>
        filter(!(get(subject) %in% in_train_set))
    )$pred
    predictions_step2 <- as.data.frame(predictions_step2)
    predictions_step2[, subject] <- not_in_train_set
    predictions_step2 <- predictions_step2 |> relocate(subject)
  }
  list(predictions_step2, not_in_train_set)
}
