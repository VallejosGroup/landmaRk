# Split a wide dataframe containing static and dynamic covariates and splits in into a dataframe with the static covariates and a list of dataframes, each associated to a dynamic covariate.

Split a wide dataframe containing static and dynamic covariates and
splits in into a dataframe with the static covariates and a list of
dataframes, each associated to a dynamic covariate.

## Usage

``` r
split_wide_df(df, ids, times, static, dynamic, measurement_name)
```

## Arguments

- df:

  A dataframe in wide format.

- ids:

  The name of the column that identifies individuals in `df`.

- times:

  The name of the column that identifies measurement times in `df`.

- static:

  A vector with the column names in `df` that store static covariates.

- dynamic:

  A vector with the column names in `df` that store dynamic covariates.

- measurement_name:

  The name for the columns where values of dynamic covariates will be
  stored.

## Value

A data frame with the static covariates, and a list of data frames, one
per dynamic covariate.
