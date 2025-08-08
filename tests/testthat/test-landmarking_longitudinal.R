initialise_longitudinal_test_ <- function(epileptic) {
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

  landmarking_object <- LandmarkAnalysis(
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
  x <- initialise_longitudinal_test_(epileptic = epileptic)
  x <- x |>
    compute_risk_sets(seq(from = 365.25, to = 5 * 365.25, by = 365.25)) |>
    fit_longitudinal(
      landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
      method = "lcmm",
      formula = value ~ treat + age + gender + learn.dis + time,
      mixture = ~ treat + age + gender + learn.dis + time,
      random = ~time,
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
      "lcmm::predictRE produced 384 predictions but expected 427 predictions.\n",
      "Probable reason: static covariates contain missing data.\n"
    )
  )
})

test_that("LOCF works as expected", {
  # Initialise LandmarkAnalysis object
  x <- initialise_longitudinal_test_(epileptic = epileptic)
  x <- x |>
    compute_risk_sets(seq(from = 365.25, to = 5 * 365.25, by = 365.25))
  expect_equal(
    {
      x <- predict_longitudinal(
        x,
        landmarks = 365.25,
        method = "locf",
        subject = "id",
        dynamic_covariates = "dose"
      )
      x@longitudinal_predictions[["365.25"]][["dose"]]
    },
    {
      locf_predictions <- x@data_dynamic[["dose"]] |>
        # Observations of dose2 before landmark time 365.25
        dplyr::filter(time <= 365.25) |>
        # Select last observation per individual
        dplyr::slice_max(time, by = id) |>
        # Join with the risk indices of individuals at risk (risk set)
        dplyr::right_join(
          data.frame(id = x@risk_sets[["365.25"]]),
          by = dplyr::join_by(id)
        ) |>
        # Select last observation
        dplyr::pull(value, name = id)
      # Mean imputation for individuals where no observations were made
      locf_predictions[is.na(locf_predictions)] <- mean(
        locf_predictions,
        na.rm = TRUE
      )
      # Sort observations according to individual id
      locf_predictions <- locf_predictions[order(as.integer(names(
        locf_predictions
      )))]
      locf_predictions
    }
  )
  expect_equal(
    {
      x <- predict_longitudinal(
        x,
        landmarks = 365.25,
        method = "locf",
        subject = "id",
        dynamic_covariates = "dose2"
      )
      x@longitudinal_predictions[["365.25"]][["dose2"]]
    },
    {
      # Query last observation per individual as in the previous test
      locf_predictions <- x@data_dynamic[["dose2"]] |>
        dplyr::filter(time <= 365.25) |>
        dplyr::slice_max(time, by = id) |>
        dplyr::right_join(
          data.frame(id = x@risk_sets[["365.25"]]),
          by = dplyr::join_by(id)
        ) |>
        dplyr::pull(value, name = id)
      # Mean imputation for individuals where no observations were made
      locf_predictions[is.na(locf_predictions)] <- names(sort(
        -table(locf_predictions)
      ))[1]
      locf_predictions <- locf_predictions[order(as.integer(names(
        locf_predictions
      )))]
      locf_predictions
    }
  )
})

test_that("longitudinal_fit raises warning for too few observations", {
  set.seed(1)
  epileptic <- epileptic |> dplyr::filter(time < 368) |> head(20)
  epileptic <- epileptic[-c(18, 19), ]

  x <- initialise_longitudinal_test_(epileptic = epileptic)

  expect_warning(
    x |>
      compute_risk_sets(seq(from = 365.25, to = 1 * 365.25, by = 365.25)) |>
      fit_longitudinal(
        landmarks = seq(from = 365.25, to = 1 * 365.25, by = 365.25),
        method = "lcmm",
        formula = value ~ treat + age + gender + learn.dis + time,
        mixture = ~ treat + age + gender + learn.dis + time,
        random = ~time,
        subject = "id",
        var.time = "time",
        ng = 2,
        dynamic_covariates = "dose"
      ),
    paste(
      "25% of the individuals have 0 or 1 observations at landmark time 365.25",
      "for longitudinal covariate dose"
    )
  )
})

test_that("predict_longitudinal works correctly with lcmm", {
  set.seed(1)

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

  x <- LandmarkAnalysis(
    data_static = static,
    data_dynamic = dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value"
  )
  x <- x |> compute_risk_sets(365.25)

  expect_warning(
    x |>
      fit_longitudinal(
        landmarks = 365.25,
        method = "lcmm",
        formula = value ~ treat + age + gender + learn.dis + time,
        mixture = ~ treat + age + gender + learn.dis + time,
        random = ~time,
        subject = "id",
        var.time = "time",
        ng = 2,
        dynamic_covariates = "dose"
      ) |>
      predict_longitudinal(
        landmarks = 365.25,
        method = "lcmm",
        subject = "id",
        var.time = "time",
        avg = FALSE,
        dynamic_covariates = "dose"
      ),
    "Individuals 28, 389, 473, have not been used in LCMM model fitting. Imputing values for those individuals"
  )
})
