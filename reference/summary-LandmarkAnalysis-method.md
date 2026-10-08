# Summarise a LandmarkAnalysis object

Summarise a LandmarkAnalysis object

## Usage

``` r
# S4 method for class 'LandmarkAnalysis'
summary(
  object,
  type = c("longitudinal", "survival"),
  landmark,
  horizon = NULL,
  dynamic_covariate = NULL
)
```

## Arguments

- object:

  An object of class
  [`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md).

- type:

  If `longitudinal`, it summarises the longitudinal submodel. If
  `survival`, it summarises the survival submodel.

- landmark:

  A numeric indicating the landmark time.

- horizon:

  For survival submodels, a numeric indicating the horizon time.

- dynamic_covariate:

  For longitudinal submodels, a character indicating the dynamic
  covariate

## Value

A summary of the desired submodel, printed to the console. Returns
`NULL` invisibly.

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
  compute_risk_sets(365.25) |>
  fit_survival(
    formula = survival::Surv(event_time, event_status) ~ treat + age,
    landmarks = 365.25,
    horizons = 2 * 365.25,
    method = "coxph"
  )
summary(x, type = "survival", landmark = 365.25, horizon = 2 * 365.25)
#> Call:
#> survival::coxph(formula = formula, data = data, model = TRUE, 
#>     x = TRUE)
#> 
#>               coef exp(coef)  se(coef)      z      p
#> treatLTG  0.111859  1.118356  0.197679  0.566 0.5715
#> age      -0.012549  0.987529  0.005414 -2.318 0.0205
#> 
#> Likelihood ratio test=5.96  on 2 df, p=0.0509
#> n= 430, number of events= 105 
# \donttest{
# Full pipeline including longitudinal submodel
x2 <- LandmarkAnalysis(
  data_static = epileptic_dfs$df_static,
  data_dynamic = epileptic_dfs$df_dynamic,
  event_indicator = "with.status",
  ids = "id", event_time = "with.time",
  times = "time", measurements = "value"
) |>
  compute_risk_sets(365.25) |>
  fit_longitudinal(
    landmarks = 365.25,
    method = "lme4",
    formula = value ~ treat + age + gender + learn.dis + (1 | id),
    dynamic_covariates = c("dose")
  ) |>
  predict_longitudinal(
    landmarks = 365.25,
    method = "lme4",
    allow.new.levels = TRUE,
    dynamic_covariates = c("dose")
  ) |>
  fit_survival(
    formula = survival::Surv(event_time, event_status) ~
      treat + age + gender + learn.dis + dose,
    landmarks = 365.25,
    horizons = 2 * 365.25,
    method = "coxph",
    dynamic_covariates = c("dose")
  ) |>
  predict_survival(landmarks = 365.25, horizons = 2 * 365.25)
summary(x2, type = "longitudinal", landmark = 365.25, dynamic_covariate = "dose")
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: value ~ treat + age + gender + learn.dis + (1 | id)
#>    Data: dataframe
#> REML criterion at convergence: 2420.099
#> Random effects:
#>  Groups   Name        Std.Dev.
#>  id       (Intercept) 0.7638  
#>  Residual             0.5132  
#> Number of obs: 1074, groups:  id, 427
#> Fixed Effects:
#>  (Intercept)      treatLTG           age       genderM  learn.disYes  
#>    2.0315626    -0.0320387    -0.0007089     0.1570440    -0.2471645  
summary(x2, type = "survival", landmark = 365.25, horizon = 2 * 365.25)
#> Call:
#> survival::coxph(formula = formula, data = data, model = TRUE, 
#>     x = TRUE)
#> 
#>                   coef exp(coef)  se(coef)      z       p
#> treatLTG      0.141183  1.151635  0.198293  0.712 0.47647
#> age          -0.015029  0.985083  0.005941 -2.530 0.01142
#> genderM      -0.022615  0.977639  0.198106 -0.114 0.90912
#> learn.disYes -0.453466  0.635422  0.436578 -1.039 0.29895
#> dose          0.375751  1.456084  0.128157  2.932 0.00337
#> 
#> Likelihood ratio test=15.45  on 5 df, p=0.008608
#> n= 430, number of events= 105 
# }
```
