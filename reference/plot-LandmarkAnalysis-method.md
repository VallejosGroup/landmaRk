# Plots longitudinal trajectories and survival curves for landmarking models.

Plots longitudinal trajectories and survival curves for landmarking
models.

## Usage

``` r
# S4 method for class 'LandmarkAnalysis'
plot(
  x,
  type = "survival",
  id = NULL,
  landmark = NULL,
  window = NULL,
  dynamic_covariate = NULL,
  avg = FALSE,
  ...
)
```

## Arguments

- x:

  An object of class
  [`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md).

- type:

  A character taking the value `'survival'` (survival curves) or
  `'longitudinal'` (model trajectories of dynamic covariates).

- id:

  The identifier for the unit (subject) whose data will be plotted.

- landmark:

  Numeric indicating a landmark time

- window:

  Numeric indicating a prediction window

- dynamic_covariate:

  A character indicating a dynamic covariate

- avg:

  A logical (by default, `FALSE`) indicating whether LCMM predictions
  are conditioned on the predicted cluster (`avg = FALSE`) or averaged
  across clusters (`avg = TRUE`). Ignored if the longitudinal model is
  not an LCMM.

- ...:

  Additional arguments passed to
  [`survminer::ggadjustedcurves()`](https://rdrr.io/pkg/survminer/man/ggadjustedcurves.html)
  for plotting survival curves.
