test_that("Error handling for fit_survival", {
  set.seed(123)

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

  landmarking_object <- Landmarking(
    data_static = static,
    data_dynamic = dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value"
  )

  x <- landmarking_object |>
    compute_risk_sets(seq(from = 365.25, to = 5 * 365.25, by = 365.25)) |>
    fit_longitudinal(
      landmarks = seq(from = 365.25, to = 2 * 365.25, by = 365.25),
      method = "lme4",
      formula = value ~ treat + age + gender + learn.dis + (time | id),
      dynamic_covariates = c("dose")
    ) |>
    predict_longitudinal(
      landmarks = seq(from = 365.25, to = 2 * 365.25, by = 365.25),
      method = "lme4",
      allow.new.levels = TRUE,
      dynamic_covariates = c("dose")
    )

  expect_error(
    x |>
      fit_survival(
        formula = Surv(event_time, event_status) ~
          treat + age + gender + learn.dis + dose,
        landmarks = seq(from = 365.25, to = 2 * 365.25, by = 365.25),
        horizons = seq(from = 2 * 365.25, to = 2 * 365.25, by = 365.25),
        method = "coxph",
        dynamic_covariates = c("dose2")
      ),
    paste(
      "Longitudinal predictions for dynamic covariate dose2 are not available at landmark time 365.25."
    )
  )

  expect_error(
    x |>
      fit_survival(
        formula = Surv(event_time, event_status) ~
          treat + age + gender + learn.dis + dose,
        landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
        horizons = seq(from = 1 * 365.25, to = 1 * 365.25, by = 365.25),
        method = "coxph",
        dynamic_covariates = c("dose")
      ),
    paste(
      "Longitudinal predictions are not available at landmark time 1095.75."
    )
  )
})
