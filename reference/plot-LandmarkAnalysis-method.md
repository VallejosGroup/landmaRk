# Plot longitudinal observations and predicted survival curve for one individual

Produces a single-panel plot with a common time axis. The left of the
landmark dashed line shows the individual's observed longitudinal
measurements; the right shows their predicted survival curve. The
summary value at the landmark that feeds into the survival sub-model is
highlighted.

## Usage

``` r
# S4 method for class 'LandmarkAnalysis'
plot(x, id, landmark, dynamic_covariate, horizon = NULL, train = TRUE, ...)
```

## Arguments

- x:

  An object of class
  [`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md).

- id:

  Identifier of the individual to plot. Must match a value in the column
  `x@ids`.

- landmark:

  Numeric landmark time.

- dynamic_covariate:

  Character name of the dynamic covariate to display.

- horizon:

  Numeric horizon time. If `NULL` (default), uses the single available
  horizon for `landmark`; errors when multiple horizons are available.

- train:

  Logical. If `TRUE` (default), uses in-sample predictions. If `FALSE`,
  uses out-of-sample predictions (requires `validation_fold > 0` in
  [`predict_survival`](https://vallejosgroup.github.io/landmaRk/reference/predict_survival.md)).

- ...:

  Additional arguments (not currently used).

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.
