# Fits the specified survival model at the landmark times and up to the horizon times specified by the user

Fits the specified survival model at the landmark times and up to the
horizon times specified by the user

## Usage

``` r
fit_survival(
  x,
  formula,
  landmarks,
  horizons,
  method,
  dynamic_covariates = c(),
  include_clusters = FALSE,
  censor_at_horizon = FALSE,
  validation_fold = 0
)
```

## Arguments

- x:

  An object of class
  [`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md).

- formula:

  A formula to be used in survival sub-model fitting.

- landmarks:

  Numeric vector of landmark times.

- horizons:

  Vector of prediction horizons up to when the survival submodel is
  fitted.

- method:

  Method for survival analysis, either "survfit" or "coxph".

- dynamic_covariates:

  Vector of time-varying covariates to be used in the survival model.

- include_clusters:

  Boolean indicating whether to propagate cluster membership to survival
  analysis.

- censor_at_horizon:

  Boolean indicating whether to censor observations at horizon times

- validation_fold:

  If positive, cross-validation fold where model is fitted. If 0
  (default), model fitting is performed on the complete dataset.

## Value

An object of class
[`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md).

## Details

### Mathematical formulation

This function estimates the conditional probability of survival to
horizon \\s+w\\, conditioned on having survived to the landmark time,
\\s\\, that is \$\$\pi_i(s+w \vert s) = P(T_i \> s+w \vert T_i \ge s,
\bar{x}\_i(s)), \$\$ where \\i\\ denotes an individual's index, \\T_i\\
is the time to event outcome for individual \\i\\ and \\\bar{x}\_i(s)\\
are the covariates observed for individual \\i\\, including the observed
history of dynamic covariates.
