test_that("auc_t uses per-subject predicted risks, not a scalar", {
  withr::local_seed(123)
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

  x <- LandmarkAnalysis(
    data_static = epileptic_dfs$df_static,
    data_dynamic = epileptic_dfs$df_dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value"
  ) |>
    compute_risk_sets(365.25) |>
    fit_survival(
      formula = survival::Surv(event_time, event_status) ~
        treat + age + gender + learn.dis,
      landmarks = 365.25,
      horizons = 2 * 365.25,
      method = "coxph"
    ) |>
    predict_survival(landmarks = 365.25, horizons = 2 * 365.25)

  result <- performance_metrics(
    x,
    landmarks = 365.25,
    horizons = 2 * 365.25,
    c_index = FALSE,
    brier = FALSE,
    auc_t = TRUE
  )

  # Manually compute the expected AUCt using per-subject predicted risks
  model_name <- "365.25-730.5"
  dataset <- x@survival_datasets[[model_name]]
  sfit <- summary(x@survival_predictions[[model_name]], times = 365.25)
  pred_risks <- 1 -
    matrix(sfit$surv, nrow = nrow(dataset), ncol = 1, byrow = TRUE)[, 1]

  score <- riskRegression::Score(
    object = list("model" = matrix(pred_risks, ncol = 1)),
    formula = survival::Surv(event_time, event_status) ~ 1,
    data = dataset,
    cause = 1,
    times = 365.25,
    cens.method = "ipcw",
    cens.model = "km"
  )$AUC$score
  expected_auct <- unname(score$AUC[score$model != "Null model"])

  expect_equal(result[, "AUCt"], expected_auct, tolerance = 1e-6)
})
