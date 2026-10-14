# Fits a survival model pooled across several landmark times (superlandmarking)

Fits a survival model pooled across several landmark times
(superlandmarking)

## Usage

``` r
fit_superlandmark(
  x,
  formula,
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

- formula:

  A formula to be used in survival submodel fitting. Must reference the
  `LM` column directly to model landmark-varying covariate effects (e.g.
  via an interaction with a function of `LM`) or a landmark-varying
  baseline hazard (via `strata(LM)`).

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

## Details

This function implements "superlandmarking": a single Cox model is
fitted to a "superdataset" subjects are present as many times as
`landmarks` they are at risk. The landmark time is available to
`formula` as the `LM` column, so that covariate effects and/or the
baseline hazard can be allowed to vary with landmark time (via an
interaction term with covariates f(LM), or `strata(LM)`).

Because each subject contributes one row per landmark at which they are
still at risk, a `cluster(@ids)` term is added to `formula`
automatically (unless already present).

## Examples

``` r
data(epileptic)
epileptic_dfs <- split_wide_df(
  epileptic,
  ids = "id", times = "time",
  static = c("with.time", "with.status", "treat", "age", "gender", "learn.dis"),
  dynamic = c("dose"),
  measurement_name = "value"
)
x <- LandmarkAnalysis(
  data_static = epileptic_dfs$df_static,
  data_dynamic = epileptic_dfs$df_dynamic,
  event_indicator = "with.status",
  ids = "id", event_time = "with.time",
  times = "time", measurements = "value"
) |>
  compute_risk_sets(c(90, 180, 365.25)) |>
  fit_superlandmark(
    formula = survival::Surv(event_time, event_status) ~
      treat * splines::ns(LM, 2) + age + strata(LM),
    landmarks = c(90, 180, 365.25),
    w = 2 * 365.25
  )
```
