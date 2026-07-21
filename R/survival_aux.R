# Check that method is coxph, finegray, survfit or a function with argument
# data, ...
.check_method_survival_predict <- function(method) {
  if (is(method)[1] == "character" && method == "survfit") {
    method <- survival::survfit
  } else if (
    is(method)[1] == "character" && method %in% c("coxph", "finegray")
  ) {
    # Leave as-is; dispatched by name in fit_survival()
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

# Warning handler shared by .fit_coxph_survival() and .fit_finegray_survival():
# raises an error instead of a warning if the underlying optimiser in
# survival::coxph() fails to converge
.coxph_convergence_handler <- function(w) {
  if (grepl("did not converge", conditionMessage(w), fixed = TRUE)) {
    stop(
      "Cox proportional hazards model failed to converge: ",
      conditionMessage(w),
      call. = FALSE
    )
  }
}

# Fits a Cox PH model via survival::coxph, raising an error instead of a
# warning if the underlying optimiser fails to converge
.fit_coxph_survival <- function(formula, data) {
  withCallingHandlers(
    survival::coxph(formula, data = data, x = TRUE, model = TRUE),
    warning = .coxph_convergence_handler
  )
}

# Fits a Fine-Gray model for the subdistribution hazard of `cause`, by
# transforming `data` with survival::finegray() and fitting a weighted Cox
# model to the result (the approach documented in ?survival::finegray).
# `data` must contain the `event_time`/`event_status` columns produced by
# .create_survival_dataframe(), with `event_status` coded as 0 = censored
# and other numeric values identifying each competing cause.
.fit_finegray_survival <- function(formula, data, cause = 1) {
  if (
    !(is.numeric(cause) &&
      length(cause) == 1L &&
      !is.na(cause) &&
      is.finite(cause))
  ) {
    stop(
      "`cause` must be a length-1, finite, non-NA numeric value.\n",
      call. = FALSE
    )
  }
  if (is.factor(data$event_status)) {
    data$event_status <- suppressWarnings(as.numeric(as.character(
      data$event_status
    )))
  }
  if (!is.numeric(data$event_status) || anyNA(data$event_status)) {
    stop(
      "`event_status` must be numeric with 0 for censoring and positive codes for causes.\n",
      call. = FALSE
    )
  }
  status_levels <- sort(unique(data$event_status))
  if (!(0 %in% status_levels)) {
    stop(
      "`event_status` must include a censoring code of 0 to fit a ",
      "Fine-Gray model.\n"
    )
  }
  if (!(cause %in% setdiff(status_levels, 0))) {
    stop(
      "`cause` (",
      cause,
      ") is not a non-censoring level of `event_status`.\n"
    )
  }

  data$event_status <- factor(
    data$event_status,
    levels = status_levels,
    labels = c("censoring", paste0("cause_", setdiff(status_levels, 0)))
  )

  finegray_data <- survival::finegray(
    formula,
    data = data,
    etype = paste0("cause_", cause)
  )

  finegray_formula <- stats::as.formula(
    paste(
      "survival::Surv(fgstart, fgstop, fgstatus) ~",
      as.character(formula)[3]
    )
  )

  # `fgwt` is looked up as a column of `finegray_data` (as in
  # ?survival::finegray's own example); it must be referenced directly rather
  # than via a local variable, since coxph() re-evaluates the `weights`
  # argument against `data`, not against this function's environment.
  withCallingHandlers(
    survival::coxph(
      finegray_formula,
      data = finegray_data,
      weights = fgwt,
      x = TRUE,
      model = TRUE
    ),
    warning = .coxph_convergence_handler
  )
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
  censor_at_horizon = FALSE,
  validation_fold = 0,
  train = TRUE
) {
  fold <- NULL
  # Create dataframe for survival analysis
  survival_df <- x@data_static[
    which(x@data_static[, x@event_time] > landmarks),
  ]
  # If censor_at_horizon=TRUE, censor observations at horizons
  if (censor_at_horizon) {
    survival_df <- survival_df |>
      mutate(
        !!sym(x@event_indicator) := ifelse(
          get(x@event_time) > horizons,
          0,
          get(x@event_indicator)
        ),
        !!sym(x@event_time) := ifelse(
          get(x@event_time) > horizons,
          horizons,
          get(x@event_time)
        )
      )
  }

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
  # Move 'baseline' to landmark time
  survival_df <- survival_df |>
    mutate(
      event_status = get(x@event_indicator),
      event_time = get(x@event_time) - landmarks
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
      if (include_clusters != FALSE && inherits(predictions, "data.frame")) {
        # Include predicted cluster membership in the training dataset and in
        # the survival formula
        survival_df <- bind_cols(
          survival_df,
          predictions
        )
        colnames(survival_df)[ncol(survival_df)] <- paste0(
          "cluster_",
          dynamic_covariate
        )
        colnames(survival_df)[ncol(survival_df) - 1] <- dynamic_covariate
      } else {
        if (is.data.frame(predictions)) {
          survival_df <- bind_cols(
            survival_df,
            predictions
          )
          colnames(survival_df)[c(
            (ncol(survival_df) - 1):ncol(survival_df)
          )] <- c(dynamic_covariate, paste0("cluster_", dynamic_covariate))
        } else {
          survival_df <- bind_cols(
            survival_df,
            matrix(
              predictions,
              ncol = 1,
              dimnames = list(rownames(survival_df), dynamic_covariate)
            )
          )
        }
      }
    }
  }

  survival_df
}
