# Helper: build tibble-based epileptic data split
.build_tibble_dfs <- function() {
  data("epileptic", package = "landmaRk", envir = environment())
  epileptic_tbl <- dplyr::as_tibble(epileptic)

  result <- split_wide_df(
    epileptic_tbl,
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

  list(
    static = dplyr::as_tibble(result$df_static),
    dynamic = list(dose = dplyr::as_tibble(result$df_dynamic[["dose"]]))
  )
}

test_that("split_wide_df accepts a tibble input", {
  data("epileptic")
  epileptic_tbl <- dplyr::as_tibble(epileptic)

  result <- split_wide_df(
    epileptic_tbl,
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

  # Outputs are data frames (tibbles satisfy is.data.frame)
  expect_true(is.data.frame(result$df_static))
  expect_true(is.data.frame(result$df_dynamic[["dose"]]))

  # Structural correctness matches the data.frame path
  df_result <- split_wide_df(
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
  expect_equal(nrow(result$df_static), nrow(df_result$df_static))
  expect_equal(
    sort(colnames(result$df_static)),
    sort(colnames(df_result$df_static))
  )
})

test_that("LandmarkAnalysis constructor accepts tibble data_static and tibble data_dynamic", {
  dfs <- .build_tibble_dfs()

  x <- LandmarkAnalysis(
    data_static = dfs$static,
    data_dynamic = dfs$dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value"
  )

  expect_s4_class(x, "LandmarkAnalysis")
  expect_true(is.data.frame(x@data_static))
  expect_true(is.data.frame(x@data_dynamic[["dose"]]))
})

test_that("Character dynamic covariate is converted to factor when data_dynamic contains tibbles", {
  dfs <- .build_tibble_dfs()

  # Make the measurement column character
  dfs$dynamic[["dose"]]$value <- as.character(dfs$dynamic[["dose"]]$value > 2)

  expect_message(
    x <- LandmarkAnalysis(
      data_static = dfs$static,
      data_dynamic = dfs$dynamic,
      event_indicator = "with.status",
      ids = "id",
      event_time = "with.time",
      times = "time",
      measurements = "value"
    ),
    "Dynamic covariate dose were coded as character. Converted to factor."
  )

  expect_equal(class(x@data_dynamic[["dose"]]$value), "factor")
})

test_that("compute_risk_sets works with tibble inputs", {
  dfs <- .build_tibble_dfs()

  x <- LandmarkAnalysis(
    data_static = dfs$static,
    data_dynamic = dfs$dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value"
  ) |>
    compute_risk_sets(landmarks = 365.25)

  expect_true(length(x@risk_sets) > 0)
  expect_true("365.25" %in% names(x@risk_sets))
})

test_that("Risk sets computed from tibble inputs match those from data.frame inputs", {
  data("epileptic")
  dfs <- .build_tibble_dfs()

  df_result <- split_wide_df(
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
  x_df <- LandmarkAnalysis(
    data_static = df_result$df_static,
    data_dynamic = df_result$df_dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value"
  ) |>
    compute_risk_sets(landmarks = 365.25)

  x_tbl <- LandmarkAnalysis(
    data_static = dfs$static,
    data_dynamic = dfs$dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value"
  ) |>
    compute_risk_sets(landmarks = 365.25)

  expect_equal(x_df@risk_sets[["365.25"]], x_tbl@risk_sets[["365.25"]])
})

test_that("compute_risk_sets drops individuals whose event_time is NA", {
  data("epileptic")
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

  # Introduce NA event_time for one individual
  na_id <- static$id[1]
  static$with.time[static$id == na_id] <- NA

  x_df <- LandmarkAnalysis(
    data_static = static,
    data_dynamic = dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value"
  ) |>
    compute_risk_sets(landmarks = 365.25)

  # The individual with NA event_time must not appear in the risk set
  expect_false(na_id %in% x_df@risk_sets[["365.25"]])

  # Same behaviour when data_static is a tibble
  static_tbl <- dplyr::as_tibble(static)
  x_tbl <- LandmarkAnalysis(
    data_static = static_tbl,
    data_dynamic = dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value"
  ) |>
    compute_risk_sets(landmarks = 365.25)

  expect_false(na_id %in% x_tbl@risk_sets[["365.25"]])
  expect_equal(x_df@risk_sets[["365.25"]], x_tbl@risk_sets[["365.25"]])
})

test_that("LOCF predict_longitudinal works with tibble inputs", {
  dfs <- .build_tibble_dfs()

  x <- LandmarkAnalysis(
    data_static = dfs$static,
    data_dynamic = dfs$dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value"
  ) |>
    compute_risk_sets(landmarks = 365.25) |>
    predict_longitudinal(
      landmarks = 365.25,
      method = "locf",
      dynamic_covariates = "dose"
    )

  predictions <- x@longitudinal_predictions[["365.25"]][["dose"]]
  expect_false(is.null(predictions))
  expect_false(any(is.na(predictions)))
})

test_that("LOCF predictions from tibble inputs match those from data.frame inputs", {
  data("epileptic")
  dfs <- .build_tibble_dfs()

  df_result <- split_wide_df(
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
  x_df <- LandmarkAnalysis(
    data_static = df_result$df_static,
    data_dynamic = df_result$df_dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value"
  ) |>
    compute_risk_sets(landmarks = 365.25) |>
    predict_longitudinal(
      landmarks = 365.25,
      method = "locf",
      dynamic_covariates = "dose"
    )

  x_tbl <- LandmarkAnalysis(
    data_static = dfs$static,
    data_dynamic = dfs$dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value"
  ) |>
    compute_risk_sets(landmarks = 365.25) |>
    predict_longitudinal(
      landmarks = 365.25,
      method = "locf",
      dynamic_covariates = "dose"
    )

  expect_equal(
    x_df@longitudinal_predictions[["365.25"]][["dose"]],
    x_tbl@longitudinal_predictions[["365.25"]][["dose"]]
  )
})

test_that("fit_longitudinal (lme4) works with tibble data_dynamic (exercises .construct_data tibble branch)", {
  dfs <- .build_tibble_dfs()

  x <- LandmarkAnalysis(
    data_static = dfs$static,
    data_dynamic = dfs$dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value",
    censor_at_landmark = TRUE
  ) |>
    compute_risk_sets(landmarks = 365.25) |>
    fit_longitudinal(
      landmarks = 365.25,
      method = "lme4",
      formula = value ~ treat + age + gender + learn.dis + (time | id),
      dynamic_covariates = "dose"
    )

  # Longitudinal model fit is available for the landmark time
  expect_false(is.null(x@longitudinal_fits[["365.25"]]))
  expect_true(is(x@longitudinal_fits[["365.25"]][["dose"]], "lmerMod"))
})
