test_that("fit_superlandmark/predict_superlandmark fit a pooled model", {
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

  landmarks <- c(365.25, 548, 730.5)
  w <- 365.25

  x <- LandmarkAnalysis(
    data_static = epileptic_dfs$df_static,
    data_dynamic = epileptic_dfs$df_dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value"
  ) |>
    compute_risk_sets(landmarks) |>
    fit_superlandmark(
      formula = Surv(event_time, event_status) ~
        treat * splines::ns(LM, 2) + age + strata(LM),
      landmarks = landmarks,
      w = w
    )

  model_name <- "365.25,548,730.5-365.25"
  expect_true(model_name %in% names(x@superlandmark_fits))
  expect_true(model_name %in% names(x@superlandmark_datasets))
  expect_s3_class(x@superlandmark_fits[[model_name]], "coxph")

  superdataset <- x@superlandmark_datasets[[model_name]]
  # One row per subject per landmark at which they were at risk.
  expect_setequal(unique(superdataset$LM), landmarks)
  # Administrative censoring at the fixed window w.
  expect_true(all(superdataset$event_time <= w))
  # cluster(id) is added automatically for robust/cluster-robust variance.
  expect_true(!is.null(x@superlandmark_fits[[model_name]]$naive.var))

  x <- predict_superlandmark(x, landmarks = landmarks, w = w)
  expect_s3_class(x@superlandmark_predictions[[model_name]], "survfit")
})

test_that("fit_superlandmark validates its inputs", {
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
    compute_risk_sets(c(365.25, 548))

  expect_error(
    fit_superlandmark(
      x,
      formula = Surv(event_time, event_status) ~ age + strata(LM),
      landmarks = 365.25,
      w = 365.25
    ),
    "at least two landmark times"
  )

  expect_error(
    fit_superlandmark(
      x,
      formula = Surv(event_time, event_status) ~ age + strata(LM),
      landmarks = c(365.25, 548),
      w = -1
    ),
    "positive numeric prediction window"
  )

  expect_error(
    fit_superlandmark(
      x,
      formula = Surv(event_time, event_status) ~ age + strata(LM),
      landmarks = c(365.25, 730.5),
      w = 365.25
    ),
    "Risk set for landmark time 730.5 has not been computed"
  )
})

test_that("predict_superlandmark errors when no model has been fitted", {
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
    compute_risk_sets(c(365.25, 548))

  expect_error(
    predict_superlandmark(x, landmarks = c(365.25, 548), w = 365.25),
    "No superlandmark model has been fitted"
  )
})

test_that(".add_cluster_term does not duplicate an existing cluster() term", {
  f <- Surv(event_time, event_status) ~ age + cluster(id) + strata(LM)
  result <- landmaRk:::.add_cluster_term(f, "id")
  expect_equal(
    length(gregexpr("cluster\\(", as.character(result)[3])[[1]]),
    1
  )
})

test_that("fit_superlandmark supports dynamic covariates from a longitudinal submodel", {
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

  landmarks <- c(365.25, 548, 730.5)
  w <- 365.25

  x <- LandmarkAnalysis(
    data_static = epileptic_dfs$df_static,
    data_dynamic = epileptic_dfs$df_dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value"
  ) |>
    compute_risk_sets(landmarks) |>
    fit_longitudinal(
      landmarks = landmarks,
      method = "lme4",
      formula = value ~ treat + age + gender + learn.dis + (time | id),
      dynamic_covariates = c("dose")
    ) |>
    predict_longitudinal(
      landmarks = landmarks,
      method = "lme4",
      allow.new.levels = TRUE,
      dynamic_covariates = c("dose")
    ) |>
    fit_superlandmark(
      formula = Surv(event_time, event_status) ~
        treat + age + dose * splines::ns(LM, 2) + strata(LM),
      landmarks = landmarks,
      w = w,
      dynamic_covariates = c("dose")
    )

  model_name <- "365.25,548,730.5-365.25"
  expect_true("dose" %in% colnames(x@superlandmark_datasets[[model_name]]))
  expect_s3_class(x@superlandmark_fits[[model_name]], "coxph")
})

test_that("plot_hazard_ratio matches directly-computed coefficients/SEs", {
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

  landmarks <- c(365.25, 548, 730.5)
  w <- 365.25

  x <- LandmarkAnalysis(
    data_static = epileptic_dfs$df_static,
    data_dynamic = epileptic_dfs$df_dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value"
  ) |>
    compute_risk_sets(landmarks) |>
    fit_survival(
      formula = Surv(event_time, event_status) ~ age,
      landmarks = landmarks,
      horizons = landmarks + w,
      method = "coxph"
    ) |>
    # No treat:LM/age:LM interaction, so the supermodel's age effect is
    # landmark-invariant by construction: a useful sanity check.
    fit_superlandmark(
      formula = Surv(event_time, event_status) ~
        treat * splines::ns(LM, 2) + age + strata(LM),
      landmarks = landmarks,
      w = w
    )

  # Crude: contrast-based log(HR) must equal the fitted coefficient directly,
  # since `age` enters that single-covariate model purely linearly.
  fit_crude <- x@survival_fits[["365.25-730.5"]]
  expect_equal(
    landmaRk:::.hazard_ratio_contrast(
      fit_crude, "age", 1,
      landmaRk:::.default_reference_row(fit_crude, "age", character(0), x@data_static)
    )$log_hr,
    unname(coef(fit_crude)["age"])
  )

  # Supermodel: age has no LM interaction, so its log(HR) must be constant
  # across landmark times and equal the fitted "age" coefficient.
  model_name <- "365.25,548,730.5-365.25"
  fit_super <- x@superlandmark_fits[[model_name]]
  hr_at_0 <- landmaRk:::.hazard_ratio_contrast(
    fit_super, "age", 1,
    data.frame(treat = factor("CBZ", levels = levels(x@data_static$treat)), LM = 0)
  )
  hr_at_5 <- landmaRk:::.hazard_ratio_contrast(
    fit_super, "age", 1,
    data.frame(treat = factor("CBZ", levels = levels(x@data_static$treat)), LM = 500)
  )
  expect_equal(hr_at_0$log_hr, hr_at_5$log_hr)
  expect_equal(hr_at_0$log_hr, unname(coef(fit_super)["age"]))

  p_both <- plot_hazard_ratio(x, "age", landmarks, w, type = "both")
  expect_s3_class(p_both, "ggplot")
  p_crude <- plot_hazard_ratio(x, "age", landmarks, w, type = "crude")
  expect_s3_class(p_crude, "ggplot")
  p_super <- plot_hazard_ratio(x, "age", landmarks, w, type = "supermodel")
  expect_s3_class(p_super, "ggplot")
})
