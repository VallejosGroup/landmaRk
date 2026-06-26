# Make predictions for time-varying covariates at specified landmark times

Make predictions for time-varying covariates at specified landmark times

## Usage

``` r
# S4 method for class 'LandmarkAnalysis'
predict_longitudinal(
  x,
  landmarks,
  method,
  dynamic_covariates,
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

  Longitudinal data analysis method used to make predictions. Either
  `"lcmm"`, `"lme4"`, `"locf"`, or a function, which can be one of two
  kinds:

  - A summary measure, like `"locf"`, computed directly from the raw
    longitudinal data and not requiring a model to have been previously
    fit with
    [`fit_longitudinal`](https://vallejosgroup.github.io/landmaRk/reference/fit_longitudinal.md).
    Such a function must have the arguments `data`, `id`, `time`,
    `value` and `landmark` (and optionally further arguments passed
    through `...`), and must return a named vector (or a two-column data
    frame) with one summary value per individual in the risk set.

  - A prediction function for a model previously fit with
    [`fit_longitudinal`](https://vallejosgroup.github.io/landmaRk/reference/fit_longitudinal.md)
    (as is the case for `"lcmm"` and `"lme4"`), where the first argument
    is the fitted model object, and which also has `newdata` and
    `subject` arguments.

- dynamic_covariates:

  Vector of time-varying covariates to be modelled as the outcome of a
  longitudinal model.

- validation_fold:

  If positive, cross-validation fold where model is fitted. If 0
  (default), model fitting is performed in the complete dataset.

- ...:

  Additional arguments passed to the prediction function (e.g. number of
  classes/clusters for lcmm).

## Value

An object of class
[`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md).
