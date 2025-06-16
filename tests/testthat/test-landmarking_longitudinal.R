test_that("LCMM works as expected", {
  set.seed(123)

  epileptic_dfs <- split_wide_df(
    epileptic,
    ids = "id",
    times = "time",
    static = c("with.time", "with.status", "treat", "age", "gender", "learn.dis"),
    dynamic = c("dose"),
    measurement_name = "value"
  )

  static <- epileptic_dfs$df_static
  dynamic <- epileptic_dfs$df_dynamic

  # dynamic$dose[sample(1:nrow(dynamic$dose), nrow(dynamic$dose) * 0.1), "value"] <- NA

  sample_missing <- sample(1:nrow(static), nrow(static) * 0.1)
  static[sample_missing, "treat"] <- NA

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
      landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
      method = "lcmm",
      formula = value ~ treat + age + gender + learn.dis,
      mixture = ~ treat + age + gender + learn.dis,
      subject = "id",
      ng = 2,
      dynamic_covariates = "dose"
    )


  expect_error(
    predict_longitudinal(x,
      landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
      method = "lcmm",
      subject = "id",
      avg = FALSE,
      dynamic_covariates = "dose"
    ),
    paste(
      "lcmm::predictY produced 387 predictions but expected 430.",
      "Probable reason: static covariates contain missing data."
    )
  )
})
