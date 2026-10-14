# Plot the hazard ratio of a covariate as a function of landmark time

Plot the hazard ratio of a covariate as a function of landmark time

## Usage

``` r
plot_hazard_ratio(
  x,
  covariate,
  landmarks,
  w,
  type = c("crude", "supermodel", "both"),
  contrast = 1,
  newdata = NULL,
  level = 0.95
)
```

## Arguments

- x:

  An object of class
  [`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md).

- covariate:

  Character name of the (numeric) covariate to compute the hazard ratio
  for.

- landmarks:

  Numeric vector of landmark times. For `type = "crude"`, a fit must
  already exist in `x@survival_fits` for each landmark (via
  [`fit_survival`](https://vallejosgroup.github.io/landmaRk/reference/fit_survival.md))
  at window `w`; for `type = "supermodel"`, a superlandmark model must
  already exist for this `landmarks`/`w` pair (via
  [`fit_superlandmark`](https://vallejosgroup.github.io/landmaRk/reference/fit_superlandmark.md)).

- w:

  Scalar prediction window, as passed to
  [`fit_survival`](https://vallejosgroup.github.io/landmaRk/reference/fit_survival.md)
  (horizons `landmarks + w`) and/or
  [`fit_superlandmark`](https://vallejosgroup.github.io/landmaRk/reference/fit_superlandmark.md).

- type:

  Whether to plot the crude per-landmark estimates, the pooled
  supermodel curve, or both.

- contrast:

  Numeric increase in `covariate` to compute the hazard ratio for
  (defaults to 1, i.e. the usual per-unit hazard ratio).

- newdata:

  A single-row `data.frame` giving values for every other covariate in
  the model (excluding `covariate` and `LM`). Defaults to the mean/first
  level of each from `x@data_static`.

- level:

  Confidence level for the pointwise intervals (default 0.95).

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Details

For every landmark time, the hazard ratio associated with a
`contrast`-unit increase in `covariate` is computed as \\\log(HR) =
(x_1 - x_0)^\top \hat\beta\\.

Any other covariate present in the formula is held fixed at the mean
(numeric) or first level (factor) of `x@data_static`, unless overridden
via `newdata`.
