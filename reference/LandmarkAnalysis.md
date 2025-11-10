# Creates an S4 class for a landmarking analysis

Creates an S4 class for a landmarking analysis

## Usage

``` r
LandmarkAnalysis(
  data_static,
  data_dynamic,
  event_indicator,
  ids,
  event_time,
  times,
  measurements,
  K = 1
)
```

## Arguments

- data_static:

  A data frame containing subject ids, static covariates,

- data_dynamic:

  Data frame in long format with subject ids, measurement values,
  measurement times and measurement name.

- event_indicator:

  Name of the column indicating event or censoring.

- ids:

  Name of the column indicating subject ids.

- event_time:

  Name of the column indicating time of the event/censoring.

- times:

  Name of the column indicating observation time in `data_dynamic`.

- measurements:

  Name of the column indicating measurement values in `data_dynamic`.

- K:

  Number of cross-validation folds (by default, 1).

## Value

An object of class `LandmarkAnalysis`
