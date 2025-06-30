initialise_longitudinal_test_ <- function() {
  set.seed(123)
  epileptic$dose2 <- as.factor(epileptic$dose > 2)

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
    dynamic = c("dose", "dose2"),
    measurement_name = "value"
  )

  static <- epileptic_dfs$df_static
  dynamic <- epileptic_dfs$df_dynamic

  sample_missing <- sample(seq_len(nrow(static)), nrow(static) * 0.1)
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



  landmarking_object
}

test_that("LCMM works as expected", {

  x <- initialise_longitudinal_test_()
  x <- x |>
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
    predict_longitudinal(
      x,
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

test_that("LOCF works as expected", {
  # Initialite Landmarking object
  x <- initialise_longitudinal_test_()
  x <- x |>
    compute_risk_sets(seq(from = 365.25, to = 5 * 365.25, by = 365.25))
  expect_equal( {
    x <- predict_longitudinal(
      x,
      landmarks = 365.25,
      method = "locf",
      subject = "id",
      dynamic_covariates = "dose"
    )
    x@longitudinal_predictions[["365.25"]][["dose"]]
  }, {
  locf_predictions <- x@data_dynamic[["dose"]] |>
    dplyr::filter(time <= 365.25) |>
    dplyr::slice_max(time, by = id)|>
    dplyr::right_join(data.frame(id = x@risk_sets[["365.25"]]), by = dplyr::join_by(id)) |>
    dplyr::pull(value, name = id)
  locf_predictions[is.na(locf_predictions)] <- mean(locf_predictions, na.rm = TRUE)
  locf_predictions <- locf_predictions[order(as.integer(names(locf_predictions)))]
  locf_predictions
  }
)
  expect_equal( {
    x <- predict_longitudinal(
      x,
      landmarks = 365.25,
      method = "locf",
      subject = "id",
      dynamic_covariates = "dose2"
    )
    x@longitudinal_predictions[["365.25"]][["dose2"]]
  }, {
    locf_predictions <- x@data_dynamic[["dose2"]] |>
      dplyr::filter(time <= 365.25) |>
      dplyr::slice_max(time, by = id)|>
      dplyr::right_join(data.frame(id = x@risk_sets[["365.25"]]), by = dplyr::join_by(id)) |>
      dplyr::pull(value, name = id)
    locf_predictions[is.na(locf_predictions)] <- names(sort(-table(locf_predictions)))[1]
    locf_predictions <- locf_predictions[order(as.integer(names(locf_predictions)))]
    locf_predictions
  }
  )
  })

