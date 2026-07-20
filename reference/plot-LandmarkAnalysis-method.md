# Plot longitudinal observations and predicted survival curve for one individual

Produces a single, self-explanatory panel with a common time axis. To
the left of the landmark dashed line, the individual's observed
longitudinal measurements are shown, together with model-based context
depending on the longitudinal sub-model used at that landmark:

- **LOCF** (or another summary measure): the last observed value is
  carried forward to the landmark (dashed segment).

- **lme4**: the population-average trajectory (fixed effects only) and
  the individual's predicted trajectory (including their predicted
  random effects) are both drawn.

- **lcmm**: the average trajectory of *every* latent cluster is drawn,
  together with the individual's own predicted trajectory (each
  cluster's fixed effects plus the individual's predicted random
  effects, averaged across clusters using the individual's posterior
  class-membership probabilities). The individual's most likely cluster
  and their posterior probability of belonging to each cluster are noted
  in the plot title/subtitle.

In every case, the value that is actually fed into the survival
sub-model is highlighted at the landmark. To the right of the landmark
dashed line, the individual's predicted survival curve is shown on a
secondary axis. A legend identifies every series.

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
