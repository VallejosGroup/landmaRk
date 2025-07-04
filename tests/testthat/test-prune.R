test_that("Prune method works", {
  # Data manipulation
  data(epileptic)

  epileptic_dfs <- split_wide_df(
    epileptic,
    ids = "id",
    times = "time",
    static = c(
      "with.time",
      "with.status",
      "treat",
      "age",
      "gender",
      "learn.dis"
    ),
    dynamic = c("dose"),
    measurement_name = "value"
  )

  static <- epileptic_dfs$df_static
  dynamic <- epileptic_dfs$df_dynamic

  landmarking_object <- LandmarkAnalysis(
    data_static = static,
    data_dynamic = dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value"
  )
  landmarking_object <- landmarking_object |>
    compute_risk_sets(landmarks = 365.25)

  landmarking_object <- landmarking_object |>
    fit_longitudinal(
      landmarks = 365.25,
      method = "lme4",
      formula = value ~ treat + age + gender + learn.dis + (time | id),
      dynamic_covariates = c("dose")
    ) |>
    predict_longitudinal(
      landmarks = 365.25,
      method = "lme4",
      allow.new.levels = TRUE,
      dynamic_covariates = c("dose")
    ) |>
    fit_survival(
      formula = Surv(event_time, event_status) ~
        treat + age + gender + learn.dis + dose,
      landmarks = 365.25,
      horizons = 2 * 365.25,
      method = "coxph",
      dynamic_covariates = c("dose")
    ) |>
    predict_survival(
      landmarks = 365.25,
      horizons = 2 * 365.25,
      method = "coxph",
      type = "survival"
    )

  expect_error(
    prune(
      landmarking_object,
      landmark = 365.26
    ),
    "The risk set at landmark time 365.26 has not been computed"
  )

  expect_length(
    unlist(list(
      landmarking_object@risk_sets,
      landmarking_object@survival_predictions,
      landmarking_object@survival_fits,
      landmarking_object@survival_datasets,
      landmarking_object@longitudinal_predictions,
      landmarking_object@longitudinal_fits
    )),
    0
  )
})
