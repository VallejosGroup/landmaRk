# Performance metrics

Computes concordance index (c-index) and Brier scores at the specified
landmark times and prediction horizons.

## Usage

``` r
# S4 method for class 'LandmarkAnalysis'
performance_metrics(
  x,
  landmarks,
  horizons,
  c_index = TRUE,
  brier = TRUE,
  auc_t = FALSE,
  train = TRUE,
  h_times = c()
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

- c_index:

  A logical. If TRUE (default), C index is reported.

- brier:

  A logical. If TRUE (default), Brier score is reported.

- auc_t:

  A logical. If TRUE, AUC_t is reported.

- train:

  A logical. If TRUE (default), performance metrics are computed in the
  training set. If FALSE, they are computed in the test set.

- h_times:

  A numeric vector of horizon times relative to the landmark times at
  which auc_t and Brier score are calculated.
