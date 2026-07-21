test_that(".fit_finegray_survival fits a convergent model without error", {
  withr::local_seed(1)
  n <- 300
  data <- data.frame(
    event_time = rexp(n, rate = 0.05),
    event_status = sample(c(0, 1, 2), n, replace = TRUE, prob = c(0.3, 0.4, 0.3)),
    age = rnorm(n, 60, 10),
    sex = factor(sample(c("M", "F"), n, replace = TRUE))
  )

  fit <- .fit_finegray_survival(
    Surv(event_time, event_status) ~ age + sex,
    data = data,
    cause = 1
  )
  expect_s3_class(fit, "coxph")
})

test_that(".fit_finegray_survival requires a censoring code of 0", {
  data <- data.frame(
    event_time = 1:10,
    event_status = rep(c(1, 2), 5),
    x = rnorm(10)
  )

  expect_error(
    .fit_finegray_survival(Surv(event_time, event_status) ~ x, data = data),
    "must include a censoring code of 0"
  )
})

test_that(".fit_finegray_survival requires a valid non-censoring `cause`", {
  data <- data.frame(
    event_time = 1:20,
    event_status = rep(c(0, 1, 2), length.out = 20),
    x = rnorm(20)
  )

  expect_error(
    .fit_finegray_survival(
      Surv(event_time, event_status) ~ x,
      data = data,
      cause = 5
    ),
    "is not a non-censoring level of `event_status`"
  )

  expect_error(
    .fit_finegray_survival(
      Surv(event_time, event_status) ~ x,
      data = data,
      cause = 0
    ),
    "is not a non-censoring level of `event_status`"
  )
})

test_that(".fit_finegray_survival raises an error when coxph fails to converge", {
  data <- data.frame(
    event_time = 1:20,
    event_status = rep(c(0, 1, 2), length.out = 20),
    x = rnorm(20)
  )

  local_mocked_bindings(
    coxph = function(...) warning("Ran out of iterations and did not converge"),
    .package = "survival"
  )

  expect_error(
    .fit_finegray_survival(Surv(event_time, event_status) ~ x, data = data),
    "Cox proportional hazards model failed to converge"
  )
})
