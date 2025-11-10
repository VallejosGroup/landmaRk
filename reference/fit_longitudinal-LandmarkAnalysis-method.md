# Fits the specified longitudinal model for time-varying covariates up to the landmark times

Fits the specified longitudinal model for time-varying covariates up to
the landmark times

## Usage

``` r
# S4 method for class 'LandmarkAnalysis'
fit_longitudinal(
  x,
  landmarks,
  method,
  formula,
  dynamic_covariates,
  validation_fold = 0,
  cores = getOption("Ncpus", 1L),
  ...
)
```

## Arguments

- x:

  An object of class
  [`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md).

- landmarks:

  A vector of Landmark times.

- method:

  Either `"lcmm"` or `"lme4"` or a function for fitting a longitudinal
  data model, where the first argument is a formula, and also has a
  `data` argument.

- formula:

  A formula to be used in longitudinal sub-model fitting.

- dynamic_covariates:

  Vector of time-varying covariates to be modelled as the outcome of a
  longitudinal model.

- validation_fold:

  If positive, cross-validation fold where model is fitted. If 0
  (default), model fitting is performed using the complete dataset.

- cores:

  Number of cores/threads to be used for parallel computation on Linux
  and MacOS. Defaults to either `options("Ncpus")` if set, or 1 (single
  threaded) otherwise. Only single-threaded computation is currently
  supported on Windows.

- ...:

  Additional arguments passed to the longitudinal model fitting function
  (e.g. number of classes/clusters for lcmm).

## Value

An object of class
[`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md).

## Details

### Parallel processing

As the longitudinal model for each landmark time is independent of the
longitudinal models for other landmark times, parallel processing can be
used to vastly speed up computation. However, due to issues with
parallel processing in R, currently only Unix-like operating systems are
supported by `landmaRk`.

## See also

[`lcmm::hlme()`](https://cecileproust-lima.github.io/lcmm/reference/hlme.html)
and [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) for
additional arguments.
