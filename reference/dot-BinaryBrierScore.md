# Binary Brier Score

Computes the Binary Brier Score (BBS) for binary outcomes at a specified
time point \\\tau\\.

## Usage

``` r
.BinaryBrierScore(predictions, time, status, tau, cause)
```

## Arguments

- predictions:

  Numeric vector of model predictions.

- time:

  Numeric vector describing the time to the event of interest or
  censoring.

- status:

  Vector of censoring status.

- tau:

  Time Brier score is evaluated.

- cause:

  Event of interest.

## Details

The BBS is a measure of the accuracy of probabilistic predictions for
binary outcomes, where the predictions are the predicted probabilities
of the event of interest occurring by time \\\tau\\. The BBS is defined
as the mean squared difference between the predicted probabilities and
the true outcome.
