# S4 class for performing a landmarking analysis

S4 class for performing a landmarking analysis

## Slots

- `landmarks`:

  A numeric vector of landmark times.

- `data_static`:

  A data frame containing subject ids, static covariates,

- `data_dynamic`:

  Data frame in long format with subject ids, measurement values,
  measurement times and measurement name.

- `event_indicator`:

  Name of the column indicating event or censoring.

- `ids`:

  Name of the column indicating subject ids.

- `times`:

  Name of the column indicating observation time in `data_dynamic`.

- `measurements`:

  Name of the column indicating measurement values in `data_dynamic`.

- `event_time`:

  Name of the column indicating time of the event/censoring.

- `risk_sets`:

  List of indices.

- `longitudinal_fits`:

  List of model fits for the specified landmark times and biomarkers.

- `longitudinal_predictions`:

  List of model predictions for the specified landmark times and
  biomarkers.

- `longitudinal_predictions_test`:

  List of out-of-sample predictions for the specified landmark times and
  biomarkers.

- `survival_datasets`:

  List of survival dataframes used in the survival submodel.

- `survival_datasets_test`:

  List of survival dataframes used for out-of-sample predictions with
  the survival submodel.

- `survival_fits`:

  List of survival model fits at each of the specified landmark times.

- `survival_predictions`:

  List of time-to-event predictions for the specified landmark times and
  prediction horizons.

- `survival_predictions_test`:

  List of out-of-sample predictions for the time-to-event outcome if K
  \> 1.

- `K`:

  Number of cross-validation folds (1 by default)

- `cv_folds`:

  Data frame associating individuals to cross-validation folds
