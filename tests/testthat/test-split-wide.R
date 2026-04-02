test_that("split_wide_df works as expected", {
  set.seed(123)

  static.correct <- c(
    "with.time",
    "with.status",
    "treat",
    "age",
    "gender",
    "learn.dis"
  )

  expect_error(
    split_wide_df(
      df = c(1, 2),
      ids = "id",
      times = "time",
      static = static.correct,
      dynamic = c("dose"),
      measurement_name = "value"
    ),
    "@df must be a data frame."
  )

  expect_error(
    split_wide_df(
      df = epileptic,
      ids = c(1, 2),
      times = "time",
      static = static.correct,
      dynamic = c("dose"),
      measurement_name = "value"
    ),
    "@ids must be a character."
  )

  expect_error(
    split_wide_df(
      df = epileptic,
      ids = c("colname1", "colname2"),
      times = "time",
      static = static.correct,
      dynamic = c("dose"),
      measurement_name = "value"
    ),
    "@ids should be of length one."
  )

  expect_error(
    split_wide_df(
      df = epileptic,
      ids = c("not.a.col"),
      times = "time",
      static = static.correct,
      dynamic = c("dose"),
      measurement_name = "value"
    ),
    "@ids must be a column name in @df."
  )

  expect_error(
    split_wide_df(
      df = epileptic,
      ids = "id",
      times = c(1, 2),
      static = static.correct,
      dynamic = c("dose"),
      measurement_name = "value"
    ),
    "@times must be a character."
  )

  expect_error(
    split_wide_df(
      df = epileptic,
      ids = "id",
      times = c("colname1", "colname2"),
      static = static.correct,
      dynamic = c("dose"),
      measurement_name = "value"
    ),
    "@times should be of length one."
  )

  expect_error(
    split_wide_df(
      df = epileptic,
      ids = "id",
      times = "time",
      static = static.correct,
      dynamic = c("dose"),
      measurement_name = 1
    ),
    "@measurement_name must be a character."
  )

  expect_error(
    split_wide_df(
      df = epileptic,
      ids = "id",
      times = "time",
      static = static.correct,
      dynamic = c("dose"),
      measurement_name = c("value", "value2")
    ),
    "@measurement_name should be of length one."
  )

  expect_error(
    split_wide_df(
      df = epileptic,
      ids = "id",
      times = "time",
      static = c(1, 2),
      dynamic = c("dose"),
      measurement_name = "value"
    ),
    "@static must be a character vector."
  )

  expect_error(
    split_wide_df(
      df = epileptic,
      ids = "id",
      times = "time",
      static = c(static.correct, "not.a.col"),
      dynamic = c("dose"),
      measurement_name = "value"
    ),
    "all elements of @static must refer to column names in @df."
  )

  expect_error(
    split_wide_df(
      df = epileptic,
      ids = "id",
      times = "time",
      static = static.correct,
      dynamic = c(1, 2),
      measurement_name = "value"
    ),
    "@dynamic must be a character vector."
  )

  expect_error(
    split_wide_df(
      df = epileptic,
      ids = "id",
      times = "time",
      static = static.correct,
      dynamic = c("dose", "not.a.col"),
      measurement_name = "value"
    ),
    "all elements of @dynamic must refer to column names in @df."
  )
})
