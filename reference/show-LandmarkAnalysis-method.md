# Displays an object of class "[`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md)"

Displays an object of class
"[`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md)"

## Usage

``` r
# S4 method for class 'LandmarkAnalysis'
show(object)
```

## Arguments

- object:

  An object of class
  [`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md).

## Value

No return value(prints a summary to the console).

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
  compute_risk_sets(365.25)
show(x)
#> Summary of LandmarkAnalysis Object:
#>   Landmarks: 365.25 
#>   Number of observations: 605 
#>   Event indicator: with.status 
#>   Event time: with.time 
#>   Risk sets: 
#>     Landmark 365.25: 430 subjects
```
