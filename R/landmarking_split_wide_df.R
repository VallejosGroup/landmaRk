#' Split a wide dataframe containing static and dynamic covariates and splits
#' in into a dataframe with the static covariates and a list of dataframes,
#' each associated to a dynamic covariate.
#'
#' @param df A dataframe in wide format.
#' @param ids The name of the column that identifies individuals in \code{df}.
#' @param times The name of the column that identifies measurement times
#'  in \code{df}.
#' @param static A vector with the column names in \code{df} that store static
#'  covariates.
#' @param dynamic A vector with the column names in \code{df} that store dynamic
#'  covariates.
#' @param measurement_name The name for the columns where values of
#' dynamic covariates will be stored.
#'
#' @returns A data frame with the static covariates, and a list of data frames,
#' one per dynamic covariate.
#' @export
#'
#' @examples
split_wide_df <- function(df, ids, times, static, dynamic, measurement_name) {
  error_str <- .build_error_str(
    df,
    ids,
    times,
    static,
    dynamic,
    measurement_name
  )
  .eval_error_str(error_str)

  df_static <- df[, c(ids, static)] |>
    unique()
  df_dynamic <- list()
  for (dyncovariate in dynamic) {
    df_dynamic[[dyncovariate]] <- df[, c(ids, times, dyncovariate)]
    colnames(df_dynamic[[dyncovariate]])[ncol(df_dynamic[[dyncovariate]])] <-
      measurement_name
  }

  return(list(
    df_static = df_static,
    df_dynamic = df_dynamic
  ))
}

.build_error_str <- function(
  df,
  ids,
  times,
  static,
  dynamic,
  measurement_name
) {
  # Build an error string for the split_wide_df function.
  error_str <- c()

  if (!is.data.frame(df)) {
    error_str <- c(error_str, "@df must be a data frame.")
  } else if (!(is(ids)[1] == "character")) {
    error_str <- c(error_str, "@ids must be a character.")
  } else if (length(ids) != 1) {
    error_str <- c(error_str, "@ids should be of length one.")
  } else if (!(is(times)[1] == "character")) {
    error_str <- c(error_str, "@times must be a character.")
  } else if (length(times) != 1) {
    error_str <- c(error_str, "@times should be of length one.")
  } else if (!(is(measurement_name)[1] == "character")) {
    error_str <- c(error_str, "@measurement_name must be a character.")
  } else if (length(measurement_name) != 1) {
    error_str <- c(error_str, "@measurement_name should be of length one.")
  } else if (!(is(static)[1] == "character")) {
    error_str <- c(error_str, "@static must be a character vector.")
  } else if (!(is(dynamic)[1] == "character")) {
    error_str <- c(error_str, "@dynamic must be a character vector.")
  } else if (!(all(static %in% colnames(df)))) {
    error_str <- c(
      error_str,
      "all elements of @static must refer to column names in @df."
    )
  } else if (!(all(dynamic %in% colnames(df)))) {
    error_str <- c(
      error_str,
      paste("all elements of @dynamic must refer to", "column names in @df.")
    )
  } else if (!(ids %in% colnames(df))) {
    error_str <- c(error_str, "@ids must be a column name in @df.")
  }
  return(error_str)
}
