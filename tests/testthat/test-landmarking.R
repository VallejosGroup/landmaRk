test_that("Validity checks for LandmarkAnalysis class work", {
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

  temp <- dynamic
  names(temp) <- NULL

  expect_error(
    LandmarkAnalysis(
      data_static = static,
      data_dynamic = temp,
      event_indicator = "with.status",
      ids = "id",
      event_time = "with.time",
      times = "time",
      measurements = "value"
    ),
    "@data_dynamic must be a named list of dataframes"
  )

  # Test: event_indicator column missing from data_static
  expect_error(
    LandmarkAnalysis(
      data_static = static,
      data_dynamic = dynamic,
      event_indicator = "missing.column",
      ids = "id",
      event_time = "with.time",
      times = "time",
      measurements = "value"
    ),
    "@event_indicator must be a column in dataframe @data_static"
  )

  # Test: ids column missing from data_static
  expect_error(
    LandmarkAnalysis(
      data_static = static,
      data_dynamic = dynamic,
      event_indicator = "with.status",
      ids = "missing.column",
      event_time = "with.time",
      times = "time",
      measurements = "value"
    ),
    "@ids must be a column in every dataframe in @data_dynamic"
  )

  # Test: event_time column missing from data_static
  expect_error(
    LandmarkAnalysis(
      data_static = static,
      data_dynamic = dynamic,
      event_indicator = "with.status",
      ids = "id",
      event_time = "missing.column",
      times = "time",
      measurements = "value"
    ),
    "@event_time must be a column in dataframe @data_static"
  )

  # Test: times column missing from data_dynamic
  expect_error(
    LandmarkAnalysis(
      data_static = static,
      data_dynamic = dynamic,
      event_indicator = "with.status",
      ids = "id",
      event_time = "with.time",
      times = "missing.column",
      measurements = "value"
    ),
    "@times must be a column in every dataframe in @data_dynamic"
  )

  # Test: measurements column missing from data_dynamic
  expect_error(
    LandmarkAnalysis(
      data_static = static,
      data_dynamic = dynamic,
      event_indicator = "with.status",
      ids = "id",
      event_time = "with.time",
      times = "time",
      measurements = "missing.column"
    ),
    "@measurements must be a column in every dataframe in @data_dynamic"
  )
})

test_that("Character covariates are converted to factor", {
  # Data manipulation
  data(epileptic)

  epileptic$gender <- as.character(epileptic$gender)
  epileptic$treat <- as.character(epileptic$treat)
  epileptic$dose <- as.character(epileptic$dose > 2)

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

  expect_message(
    x <- LandmarkAnalysis(
      data_static = static,
      data_dynamic = dynamic,
      event_indicator = "with.status",
      ids = "id",
      event_time = "with.time",
      times = "time",
      measurements = "value"
    ),
    "Static covariates treat, gender were coded as characters. Converted to factors."
  )

  expect_equal(class(x@data_static$treat), "factor")
  expect_equal(class(x@data_static$gender), "factor")

  dynamic[["dose"]]$value <- as.character(dynamic[["dose"]]$value > 2)

  expect_message(
    x <- LandmarkAnalysis(
      data_static = static,
      data_dynamic = dynamic,
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

test_that("Prune individuals from risk set works", {
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

  x <- LandmarkAnalysis(
    data_static = static,
    data_dynamic = dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value"
  ) |>
    compute_risk_sets(
      landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25)
    )

  expect_message(
    x |> prune_risk_sets(365.25, c(600, 601)),
    "Removing 2 individuals from risk set for landmark time 365.25"
  )

  expect_warning(
    x |> prune_risk_sets(365.25, c(603, 604, 605)),
    "A total of  3 in @individuals are not in the risk set for landmark time 365.25"
  )

  expect_equal(
    suppressWarnings(prune_risk_sets(x, 365.25, c(603, 604, 605))@risk_sets[[
      "365.25"
    ]]),
    setdiff(x@risk_sets[["365.25"]], c(603, 604, 605))
  )
})

test_that("Warning is raised with few observations per individual", {
  # Data manipulation
  data(epileptic)

  epileptic <- epileptic |> filter(id < 10)
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

  expect_warning(
    x |>
      compute_risk_sets(
        landmarks = 365.25,
        .warn_when_less_than = 2
      ),
    "The following individuals have less than 2 observations recorded prior to landmark time 365.25 for dynamic covariate dose: 7"
  )
})
