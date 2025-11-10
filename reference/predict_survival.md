# Make predictions for time-to-event outcomes at specified horizon times

Make predictions for time-to-event outcomes at specified horizon times

## Usage

``` r
predict_survival(
  x,
  landmarks,
  horizons,
  method,
  dynamic_covariates = c(),
  include_clusters = FALSE,
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

- horizons:

  Vector of prediction horizons up to when the survival submodel is
  fitted.

- method:

  R function that is used to make predictions

- dynamic_covariates:

  Vector of time-varying covariates to be used in the survival model.

- include_clusters:

  Boolean indicating whether to propagate cluster membership to survival
  analysis.

- validation_fold:

  If positive, cross-validation fold where model is fitted. If 0
  (default), model fitting is performed on the complete dataset.

- ...:

  Additional arguments passed to the prediction function (e.g. number of
  classes/clusters for lcmm).

## Value

An object of class
[`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md).
