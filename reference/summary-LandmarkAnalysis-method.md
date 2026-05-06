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

A summary of the desired submodel

## Examples

``` r
# \donttest{
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
  fit_longitudinal(
    landmarks = 365.25,
    method = "lme4",
    formula = value ~ treat + age + gender + learn.dis + (time | id),
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
summary(x, type = "longitudinal", landmark = 365.25, dynamic_covariate = "dose")
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: value ~ treat + age + gender + learn.dis + (time | id)
#>    Data: dataframe
#> REML criterion at convergence: 2246.377
#> Random effects:
#>  Groups   Name        Std.Dev. Corr  
#>  id       (Intercept) 0.713682       
#>           time        0.003222 -0.22 
#>  Residual             0.358703       
#> Number of obs: 1074, groups:  id, 427
#> Fixed Effects:
#>  (Intercept)      treatLTG           age       genderM  learn.disYes  
#>    1.9585547    -0.1244086    -0.0006828     0.1257524    -0.2773500  
#> optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 2 lme4 warnings 
summary(x, type = "survival", landmark = 365.25, horizon = 2 * 365.25)
#> Call:
#> survival::coxph(formula = formula, data = x@survival_datasets[[paste0(landmarks, 
#>     "-", horizons)]], model = TRUE, x = TRUE)
#> 
#>                   coef exp(coef)  se(coef)      z        p
#> treatLTG      0.110243  1.116549  0.197993  0.557 0.577664
#> age          -0.014894  0.985216  0.005952 -2.502 0.012341
#> genderM      -0.022325  0.977923  0.198184 -0.113 0.910311
#> learn.disYes -0.488076  0.613806  0.436477 -1.118 0.263475
#> dose          0.276647  1.318700  0.082287  3.362 0.000774
#> 
#> Likelihood ratio test=17.86  on 5 df, p=0.003126
#> n= 430, number of events= 105 
# }
```
