# Check that method is coxph, survfit or a function with argument data, ...
.check_method_survival_predict <- function(method) {
  if (is(method)[1] == "character" && method == "survfit") {
    method <- survival::survfit
  } else if (is(method)[1] == "character" && method == "coxph") {
    method <- "coxph"
  } else if (!(is(method)[1] == "function")) {
    stop(
      "Argument ",
      method,
      " must be a function",
      "\n"
    )
  }
  if (
    (is(method)[1] == "function") &&
      !("data" %in% names(as.list(args(method))))
  ) {
    stop(
      "Argument ",
      method,
      " must be a function, and data must be an argument to that function",
      "\n"
    )
  }
  method
}

# Check risk set is available
.check_riskset_survival <- function(x, landmarks) {
  if (!(landmarks %in% x@landmarks)) {
    stop(
      "Risk set for landmark time ",
      landmarks,
      " has not been computed\n"
    )
  }
}

.check_predictions_available_survival <- function(
  x,
  landmarks,
  dynamic_covariates
) {
  landmarks <- as.character(landmarks)
  if (!(landmarks %in% names(x@longitudinal_predictions))) {
    stop(
      "Longitudinal predictions are not available at landmark time ",
      landmarks,
      "."
    )
  }

  if (length(dynamic_covariates) > 0) {
    for (dynamic_covariate in dynamic_covariates) {
      if (
        !(dynamic_covariate %in% names(x@longitudinal_predictions[[landmarks]]))
      ) {
        stop(
          "Longitudinal predictions for dynamic covariate ",
          dynamic_covariate,
          " are not available at landmark time ",
          landmarks,
          "."
        )
      }
    }
  }
}

.create_survival_dataframe <- function(
  x,
  landmarks,
  horizons,
  dynamic_covariates = c(),
  include_clusters = FALSE,
  validation_fold = 0,
  train = TRUE
) {
  fold <- NULL
  # Create dataframe for survival analysis
  survival_df <- x@survival_datasets[[paste0(
    landmarks,
    "-",
    horizons
  )]] <- x@data_static[
    which(x@data_static[, x@event_time] >= landmarks),
  ]
  # Select individuals for train or for test
  if (train) {
    survival_df <- survival_df |>
      inner_join(
        x@cv_folds |>
          filter(fold != validation_fold) |>
          select(x@ids),
        by = x@ids
      )
  } else {
    survival_df <- survival_df |>
      inner_join(
        x@cv_folds |>
          filter(fold == validation_fold) |>
          select(x@ids),
        by = x@ids
      )
  }
  # Censor observations past the horizon time
  survival_df <- survival_df |>
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
    # Add predicted values for dynamic covariates to survival training dataset
    for (dynamic_covariate in dynamic_covariates) {
      if (train) {
        predictions <- x@longitudinal_predictions[[as.character(landmarks)]][[
          dynamic_covariate
        ]]
      } else {
        predictions <- x@longitudinal_predictions_test[[as.character(
          landmarks
        )]][[
          dynamic_covariate
        ]]
      }
      survival_df <- bind_cols(
        survival_df,
        predictions
      )
      if (include_clusters && inherits(predictions, "data.frame")) {
        # Include predicted cluster membership in the training dataset and in the survival formula
        colnames(survival_df)[ncol(survival_df)] <- paste0(
          "cluster_",
          dynamic_covariate
        )
        colnames(survival_df)[ncol(survival_df) - 1] <- dynamic_covariate
      } else {
        colnames(survival_df)[ncol(survival_df)] <- dynamic_covariate
      }
    }
  }

  survival_df
}
