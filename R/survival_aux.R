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
