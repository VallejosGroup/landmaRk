test_that("Validity checks for Landmarking class work", {
  # Data manipulation
  data(epileptic)

  epilectic_dfs <- split_wide_df(
    epileptic,
    ids = "id",
    times = "time",
    static = c("with.time", "with.status", "treat", "age", "gender", "learn.dis"),
    dynamic = c("dose"),
    measurement_name = "value"
  )

  static <- epilectic_dfs$df_static
  dynamic <- epilectic_dfs$df_dynamic

  temp <- dynamic
  names(temp) <- NULL

  expect_error(
    Landmarking(
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
    Landmarking(
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
    Landmarking(
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
    Landmarking(
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
    Landmarking(
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
    Landmarking(
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
