test_that("Hist is available without attaching prodlim", {
  # Regression test for issue #149
  # Ensures that Hist from prodlim can be used in formulas without
  # requiring users to call library(prodlim)
  
  # Verify prodlim is NOT attached (not in search path)
  expect_false("package:prodlim" %in% search())
  
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
  
  # Create LandmarkAnalysis object and run through pipeline
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
    compute_risk_sets(landmarks = 365.25) |>
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
  
  # Test performance_metrics - this calls riskRegression::Score
  # which evaluates formulas that may contain Hist()
  # Without the import, this would fail with "could not find function 'Hist'"
  expect_no_error(
    metrics <- performance_metrics(
      landmarking_object,
      landmarks = 365.25,
      horizons = 2 * 365.25,
      c_index = TRUE,
      brier = TRUE,
      train = TRUE
    )
  )
  
  expect_true(is.data.frame(metrics) || is.matrix(metrics))
  expect_true("cindex" %in% colnames(metrics))
  expect_true("Brier" %in% colnames(metrics))
  
  # Verify prodlim is still NOT attached after the test
  expect_false("package:prodlim" %in% search())
})
