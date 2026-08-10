test_that(".fit_coxph_survival fits a convergent model without error", {
  withr::local_seed(1)
  n <- 200
  data <- data.frame(
    time = rexp(n),
    status = rbinom(n, 1, 0.7),
    x = rnorm(n)
  )

  fit <- landmaRk:::.fit_coxph_survival(Surv(time, status) ~ x, data = data)
  expect_s3_class(fit, "coxph")
})

test_that(".fit_coxph_survival raises an error when coxph fails to converge", {
  data <- data.frame(time = 1:10, status = rep(1, 10), x = rnorm(10))

  local_mocked_bindings(
    coxph = function(...) warning("Ran out of iterations and did not converge"),
    .package = "survival"
  )

  expect_error(
    landmaRk:::.fit_coxph_survival(Surv(time, status) ~ x, data = data),
    "Cox proportional hazards model failed to converge"
  )
})

test_that(".fit_coxph_survival propagates unrelated warnings without erroring", {
  data <- data.frame(time = 1:10, status = rep(1, 10), x = rnorm(10))

  local_mocked_bindings(
    coxph = function(...) {
      warning("some unrelated warning")
      "fit"
    },
    .package = "survival"
  )

  expect_warning(
    result <- landmaRk:::.fit_coxph_survival(
      Surv(time, status) ~ x,
      data = data
    ),
    "some unrelated warning"
  )
  expect_equal(result, "fit")
})

test_that("fit_survival raises an error when the Cox PH model fails to converge", {
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

  landmarking_object <- LandmarkAnalysis(
    data_static = epileptic_dfs$df_static,
    data_dynamic = epileptic_dfs$df_dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value"
  ) |>
    compute_risk_sets(landmarks = 365.25)

  local_mocked_bindings(
    coxph = function(...) warning("Ran out of iterations and did not converge"),
    .package = "survival"
  )

  expect_error(
    landmarking_object |>
      fit_survival(
        formula = Surv(event_time, event_status) ~
          treat + age + gender + learn.dis,
        landmarks = 365.25,
        horizons = 2 * 365.25,
        method = "coxph"
      ),
    "Cox proportional hazards model failed to converge"
  )
})
