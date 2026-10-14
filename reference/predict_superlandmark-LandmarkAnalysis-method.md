# Make predictions from a superlandmark survival model

Make predictions from a superlandmark survival model

## Usage

``` r
# S4 method for class 'LandmarkAnalysis'
predict_superlandmark(
  x,
  landmarks,
  w,
  dynamic_covariates = c(),
  include_clusters = FALSE,
  validation_fold = 0
)
```

## Arguments

- x:

  An object of class
  [`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md).

- landmarks:

  Numeric vector of at least two landmark times, for which risk sets
  must already have been computed via
  [`compute_risk_sets`](https://vallejosgroup.github.io/landmaRk/reference/compute_risk_sets.md).

- w:

  Scalar prediction window. Each landmark time \\s\\ in `landmarks` uses
  the (fixed) horizon \\s+w\\.

- dynamic_covariates:

  Vector of time-varying covariates to be used in the survival model.

- include_clusters:

  Boolean indicating whether to propagate cluster membership to survival
  analysis.

- validation_fold:

  If positive, cross-validation fold where model is fitted. If 0
  (default), model fitting is performed on the complete dataset.

## Value

An object of class
[`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md).
