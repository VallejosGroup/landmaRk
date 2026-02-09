# Cross-validation with the landmaRk package

## Setup

In addition to the `landmaRk` package, we will also use `tidyverse`.

``` r
set.seed(123)
library(landmaRk)
library(tidyverse)
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.2.0     ✔ readr     2.1.6
#> ✔ forcats   1.0.1     ✔ stringr   1.6.0
#> ✔ ggplot2   4.0.2     ✔ tibble    3.3.1
#> ✔ lubridate 1.9.5     ✔ tidyr     1.3.2
#> ✔ purrr     1.2.1     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
library(survival)
library(prodlim)
```

## Example: `aids` data

As in the first vignette, we use the epileptic dataset to perform
landmarking analysis of time-to-event data with time-varying covariates.
Here is the structure of the dataset.

``` r
library(JMbayes2)
#> Loading required package: nlme
#> 
#> Attaching package: 'nlme'
#> The following object is masked from 'package:dplyr':
#> 
#>     collapse
#> Loading required package: GLMMadaptive
#> Loading required package: splines
#> 
#> Attaching package: 'JMbayes2'
#> The following object is masked from 'package:GLMMadaptive':
#> 
#>     mixed_model
#> The following object is masked from 'package:nlme':
#> 
#>     lme
#> The following object is masked from 'package:survival':
#> 
#>     coxph
data("aids")
aids$patient <- as.numeric(aids$patient)
str(aids)
#> 'data.frame':    1405 obs. of  12 variables:
#>  $ patient: num  1 1 1 2 2 2 2 3 3 3 ...
#>  $ Time   : num  17 17 17 19 19 ...
#>  $ death  : int  0 0 0 0 0 0 0 1 1 1 ...
#>  $ CD4    : num  10.68 8.43 9.43 6.32 8.12 ...
#>  $ obstime: int  0 6 12 0 6 12 18 0 2 6 ...
#>  $ drug   : Factor w/ 2 levels "ddC","ddI": 1 1 1 2 2 2 2 2 2 2 ...
#>  $ gender : Factor w/ 2 levels "female","male": 2 2 2 2 2 2 2 1 1 1 ...
#>  $ prevOI : Factor w/ 2 levels "noAIDS","AIDS": 2 2 2 1 1 1 1 2 2 2 ...
#>  $ AZT    : Factor w/ 2 levels "intolerance",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ start  : int  0 6 12 0 6 12 18 0 2 6 ...
#>  $ stop   : num  6 12 17 6 12 ...
#>  $ event  : num  0 0 0 0 0 0 0 0 0 1 ...
```

## Initialising the landmarking analysis for cross-validation

As in the introductory vignette, the dataset `epileptic` comes in wide
format. We split it into two dataframes, one for static covariates and
one for dynamic covariates.

``` r
# DF with Static covariates
aids_dfs <- split_wide_df(
  aids,
  ids = "patient",
  times = "obstime",
  static = c("Time", "death", "drug", "gender", "prevOI"),
  dynamic = c("CD4"),
  measurement_name = "value"
)
static <- aids_dfs$df_static
head(static)
#>    patient  Time death drug gender prevOI
#> 1        1 16.97     0  ddC   male   AIDS
#> 4        2 19.00     0  ddI   male noAIDS
#> 8        3 18.53     1  ddI female   AIDS
#> 11       4 12.70     0  ddC   male   AIDS
#> 15       5 15.13     0  ddI   male   AIDS
#> 19       6  1.90     1  ddC female   AIDS
```

``` r
# DF with Dynamic covariates
dynamic <- aids_dfs$df_dynamic
head(dynamic[["CD4"]])
#>   patient obstime     value
#> 1       1       0 10.677078
#> 2       1       6  8.426150
#> 3       1      12  9.433981
#> 4       2       0  6.324555
#> 5       2       6  8.124038
#> 6       2      12  4.582576
```

As in the introductory vignette, we create an object of class
`LandmarkAnalysis`, using the helper function of the same name. We now
use the optional argument `K` to specify the number of cross-validations
folds. In this example, we request five cross validation folds.

We then calculate the risk sets using
[`compute_risk_sets()`](https://vallejosgroup.github.io/landmaRk/reference/compute_risk_sets.md).

``` r
landmarking_object <- LandmarkAnalysis(
  data_static = static,
  data_dynamic = dynamic,
  event_indicator = "death",
  ids = "patient",
  event_time = "Time",
  times = "obstime",
  measurements = "value",
  K = 5
)

landmarking_object <- landmarking_object |>
  compute_risk_sets(landmarks = c(6, 8))
```

## Performance evaluation in hold-out dataset

Now that we have split the dataset into `K=5` parts for
cross-validation, we can use one of the five parts as test set and the
remaining four parts as the training set. To do so, use the argument
`validation_fold` to indicate the that you want to use as test set when
calling
[`fit_longitudinal()`](https://vallejosgroup.github.io/landmaRk/reference/fit_longitudinal.md),
[`predict_longitudinal()`](https://vallejosgroup.github.io/landmaRk/reference/predict_longitudinal.md),
[`fit_survival()`](https://vallejosgroup.github.io/landmaRk/reference/fit_survival.md)
and
[`predict_survival()`](https://vallejosgroup.github.io/landmaRk/reference/predict_survival.md).

``` r
landmarking_object <- landmarking_object |>
  fit_longitudinal(
    landmarks = c(6, 8),
    method = "lme4",
    formula = value ~ prevOI + obstime + (obstime | patient),
    dynamic_covariates = c("CD4"),
    validation_fold = 5
  ) |>
  predict_longitudinal(
    landmarks = c(6, 8),
    method = "lme4",
    dynamic_covariates = c("CD4"),
    validation_fold = 5,
    allow.new.levels = TRUE
  ) |>
  fit_survival(
    formula = Surv(event_time, event_status) ~ drug,
    landmarks = c(6, 8),
    horizons = 12 + c(6, 8),
    method = "coxph",
    dynamic_covariates = c("CD4"),
    validation_fold = 5
  ) |>
  predict_survival(
    landmarks = c(6, 8),
    horizons = 12 + c(6, 8),
    method = "coxph",
    type = "lp",
    validation_fold = 5
  )
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> New names:
#> New names:
#> • `` -> `...9`
```

We can also use [`summary()`](https://rdrr.io/r/base/summary.html) to
print parameter estimates. Note that this time the sample size is
smaller, because we have held out part of the original dataset for
validation.

``` r
summary(landmarking_object, type = "longitudinal", landmark = 6, dynamic_covariate = "CD4")
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: value ~ prevOI + obstime + (obstime | patient)
#>    Data: dataframe
#> REML criterion at convergence: 4263.201
#> Random effects:
#>  Groups   Name        Std.Dev. Corr 
#>  patient  (Intercept) 4.1584        
#>           obstime     0.2392   -0.09
#>  Residual             1.6223        
#> Number of obs: 853, groups:  patient, 321
#> Fixed Effects:
#> (Intercept)   prevOIAIDS      obstime  
#>     10.3201      -4.3287      -0.1757
```

``` r
summary(landmarking_object, type = "survival", landmark = 6, horizon = 18)
#> Call:
#> survival::coxph(formula = formula, data = x@survival_datasets[[paste0(landmarks, 
#>     "-", horizons)]], model = TRUE, x = TRUE)
#> 
#>           coef exp(coef) se(coef)     z   p
#> drugddI 0.2081    1.2313   0.2007 1.037 0.3
#> 
#> Likelihood ratio test=1.08  on 1 df, p=0.299
#> n= 321, number of events= 100
```

Here are the in-sample performance metrics:

``` r
performance_metrics(
  landmarking_object,
  landmarks = c(6, 8),
  horizons = c(18, 20), 
  auc_t = TRUE, c_index = FALSE,
  h_times = c(3, 6, 12)
)
#>      landmark horizon   Brier(9) Brier(12) Brier(18)    AUC(3)    AUC(6)
#> 6-18        6      18 0.07967156 0.1630371 0.2313983 0.5194303 0.5381293
#> 8-20        8      20 0.10531988 0.1699368 0.2449021 0.5118598 0.5353575
#>        AUC(12)
#> 6-18 0.4853891
#> 8-20 0.4644064
```

Out-of-sample performance metrics can be obtained by specifying
`train = FALSE`:

``` r
performance_metrics(
  landmarking_object,
  landmarks = c(6, 8),
  horizons = c(18, 20), 
  auc_t = TRUE, c_index = FALSE,
  h_times = c(3, 6, 12),
  train = FALSE
)
#>      landmark horizon   Brier(9) Brier(12) Brier(18)    AUC(3)    AUC(6)
#> 6-18        6      18 0.07655582 0.1744548 0.2385842 0.6400376 0.6126645
#> 8-20        8      20 0.08189975 0.1732424 0.2365391 0.6458753 0.6614730
#>        AUC(12)
#> 6-18 0.5646015
#> 8-20 0.4347708
```

## Cross-validation

Now, we can embed the above pipeline in a for loop to perform
cross-validation:

``` r
landmarking_object <- LandmarkAnalysis(
  data_static = static,
  data_dynamic = dynamic,
  event_indicator = "death",
  ids = "patient",
  event_time = "Time",
  times = "obstime",
  measurements = "value",
  K = 5
)

landmarking_object <- landmarking_object |>
  compute_risk_sets(landmarks = c(6, 8))
```

``` r
metrics <- list()
for (k in 1:5) {
  metrics[[k]] <- landmarking_object |>
    fit_longitudinal(
    landmarks = c(6, 8),
    method = "lme4",
    formula = value ~ prevOI + obstime + (obstime | patient),
    dynamic_covariates = c("CD4"),
    validation_fold = k
  ) |>
  predict_longitudinal(
    landmarks = c(6, 8),
    method = "lme4",
    dynamic_covariates = c("CD4"),
    validation_fold = k,
    allow.new.levels = TRUE
  ) |>
  fit_survival(
    formula = Surv(event_time, event_status) ~ drug,
    landmarks = c(6, 8),
    horizons = 12 + c(6, 8),
    method = "coxph",
    dynamic_covariates = c("CD4"),
    validation_fold = k
  ) |>
  predict_survival(
    landmarks = c(6, 8),
    horizons = 12 + c(6, 8),
    method = "coxph",
    type = "lp",
    validation_fold = k
  ) |>
    performance_metrics(
      landmarks = c(6, 8),
      horizons = c(18, 20), 
      auc_t = TRUE, brier = TRUE, c_index = FALSE,
      h_times = c(3, 6, 12)
    )
}
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> New names:
#> New names:
#> • `` -> `...9`
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> New names:
#> New names:
#> • `` -> `...9`
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> New names:
#> New names:
#> • `` -> `...9`
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> New names:
#> New names:
#> • `` -> `...9`
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> New names:
#> New names:
#> • `` -> `...9`

metrics
#> [[1]]
#>      landmark horizon   Brier(9) Brier(12) Brier(18)    AUC(3)    AUC(6)
#> 6-18        6      18 0.07430174 0.1628066 0.2331222 0.5165677 0.5332049
#> 8-20        8      20 0.10693845 0.1686616 0.2472218 0.5129526 0.5464116
#>        AUC(12)
#> 6-18 0.4909934
#> 8-20 0.4711122
#> 
#> [[2]]
#>      landmark horizon   Brier(9) Brier(12) Brier(18)    AUC(3)    AUC(6)
#> 6-18        6      18 0.07878887 0.1649606 0.2357098 0.5245006 0.5430452
#> 8-20        8      20 0.10599282 0.1717609 0.2397350 0.5231481 0.5536476
#>        AUC(12)
#> 6-18 0.5008851
#> 8-20 0.4908798
#> 
#> [[3]]
#>      landmark horizon   Brier(9) Brier(12) Brier(18)    AUC(3)    AUC(6)
#> 6-18        6      18 0.07071204 0.1598846 0.2288635 0.5448391 0.5714684
#> 8-20        8      20 0.09567283 0.1690807 0.2444550 0.5446251 0.5772778
#>        AUC(12)
#> 6-18 0.5021795
#> 8-20 0.4210204
#> 
#> [[4]]
#>      landmark horizon   Brier(9) Brier(12) Brier(18)    AUC(3)    AUC(6)
#> 6-18        6      18 0.08955414 0.1686272 0.2363873 0.5502151 0.5302593
#> 8-20        8      20 0.09355146 0.1656261 0.2322180 0.5037594 0.5340163
#>        AUC(12)
#> 6-18 0.4867656
#> 8-20 0.4807579
#> 
#> [[5]]
#>      landmark horizon   Brier(9) Brier(12) Brier(18)    AUC(3)    AUC(6)
#> 6-18        6      18 0.08130580 0.1681432 0.2359809 0.5782629 0.5907356
#> 8-20        8      20 0.09911599 0.1739465 0.2599817 0.5934627 0.5997129
#>        AUC(12)
#> 6-18 0.5300017
#> 8-20 0.4345596
```
