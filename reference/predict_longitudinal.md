# Make predictions for time-varying covariates at specified landmark times

Make predictions for time-varying covariates at specified landmark times

## Usage

``` r
predict_longitudinal(
  x,
  landmarks,
  method,
  dynamic_covariates,
  censor_at_horizon = FALSE,
  validation_fold = 0,
  ...
)
```

## Arguments

- x:

  An object of class
  [`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md).

- landmarks:

  A numeric vector of landmark times.

- method:

  Longitudinal data analysis method used to make predictions

- dynamic_covariates:

  Vector of time-varying covariates to be modelled as the outcome of a
  longitudinal model.

- censor_at_horizon:

  Boolean indicating whether to censor observations at horizon times

- validation_fold:

  If positive, cross-validation fold where model is fitted. If 0
  (default), model fitting is performed in the complete dataset.

- ...:

  Additional arguments passed to the prediction function (e.g. number of
  classes/clusters for lcmm).

## Value

An object of class
[`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md).
