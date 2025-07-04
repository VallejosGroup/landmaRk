% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/landmarking_split_wide_df.R
\name{split_wide_df}
\alias{split_wide_df}
\title{Split a wide dataframe containing static and dynamic covariates and splits
in into a dataframe with the static covariates and a list of dataframes,
each associated to a dynamic covariate.}
\usage{
split_wide_df(df, ids, times, static, dynamic, measurement_name)
}
\arguments{
\item{df}{A dataframe in wide format.}

\item{ids}{The name of the column that identifies individuals in \code{df}.}

\item{times}{The name of the column that identifies measurement times
in \code{df}.}

\item{static}{A vector with the column names in \code{df} that store static
covariates.}

\item{dynamic}{A vector with the column names in \code{df} that store dynamic
covariates.}

\item{measurement_name}{The name for the columns where values of
dynamic covariates will be stored.}
}
\value{
A data frame with the static covariates, and a list of data frames,
one per dynamic covariate.
}
\description{
Split a wide dataframe containing static and dynamic covariates and splits
in into a dataframe with the static covariates and a list of dataframes,
each associated to a dynamic covariate.
}
