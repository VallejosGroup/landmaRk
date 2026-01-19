# Cross-validation with the landmaRk package

## Setup

In addition to the `landmaRk` package, we will also use `tidyverse`.

``` r
set.seed(123)
library(landmaRk)
library(tidyverse)
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.4     ✔ readr     2.1.6
#> ✔ forcats   1.0.1     ✔ stringr   1.6.0
#> ✔ ggplot2   4.0.1     ✔ tibble    3.3.1
#> ✔ lubridate 1.9.4     ✔ tidyr     1.3.2
#> ✔ purrr     1.2.1     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
library(survival)
library(prodlim)
```

## Example: epileptic data

As in the first vignette, we use the epileptic dataset to perform
landmarking analysis of time-to-event data with time-varying covariates.
Here is the structure of the dataset.

``` r
data("epileptic")
str(epileptic)
#> 'data.frame':    2797 obs. of  9 variables:
#>  $ id         : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ time       : int  86 119 268 451 535 770 1515 1829 2022 2194 ...
#>  $ with.time  : int  2400 2400 2400 2400 2400 2400 2400 2400 2400 2400 ...
#>  $ with.status: int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ dose       : num  2 2 2 2.67 2.67 2.67 2.67 2.67 3.33 2.67 ...
#>  $ treat      : Factor w/ 2 levels "CBZ","LTG": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ age        : num  75.7 75.7 75.7 75.7 75.7 ...
#>  $ gender     : Factor w/ 2 levels "F","M": 2 2 2 2 2 2 2 2 2 2 ...
#>  $ learn.dis  : Factor w/ 2 levels "No","Yes": 1 1 1 1 1 1 1 1 1 1 ...
```

## Initialising the landmarking analysis for cross-validation

As in the introductory vignette, the dataset `epileptic` comes in wide
format. We split it into two dataframes, one for static covariates and
one for dynamic covariates.

``` r
# DF with Static covariates
epileptic_dfs <- split_wide_df(
  epileptic,
  ids = "id",
  times = "time",
  static = c("with.time", "with.status", "treat", "age", "gender", "learn.dis"),
  dynamic = c("dose"),
  measurement_name = "value"
)
static <- epileptic_dfs$df_static
head(static)
#>    id with.time with.status treat   age gender learn.dis
#> 1   1      2400           0   CBZ 75.67      M        No
#> 12  2      2324           0   LTG 32.96      M        No
#> 25  3       802           0   LTG 29.31      M        No
#> 29  4      2364           0   CBZ 44.59      M        No
#> 42  5       821           1   LTG 40.61      F        No
#> 45  6      2237           0   LTG 28.06      M       Yes
```

``` r
# DF with Dynamic covariates
dynamic <- epileptic_dfs$df_dynamic
head(dynamic[["dose"]])
#>   id time value
#> 1  1   86  2.00
#> 2  1  119  2.00
#> 3  1  268  2.00
#> 4  1  451  2.67
#> 5  1  535  2.67
#> 6  1  770  2.67
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
  event_indicator = "with.status",
  ids = "id",
  event_time = "with.time",
  times = "time",
  measurements = "value",
  K = 5
)

landmarking_object <- landmarking_object |>
  compute_risk_sets(
    landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25)
  )
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
    landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
    method = "lme4",
    formula = value ~ treat + age + gender + learn.dis + (time | id),
    dynamic_covariates = c("dose"),
    validation_fold = 5
  ) |>
  predict_longitudinal(
    landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
    method = "lme4",
    dynamic_covariates = c("dose"),
    validation_fold = 5,
    allow.new.levels = TRUE
  ) |>
  fit_survival(
    formula = Surv(event_time, event_status) ~
      treat + age + gender + learn.dis + dose,
    landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
    horizons = seq(from = 2 * 365.25, to = 6 * 365.25, by = 365.25),
    method = "coxph",
    dynamic_covariates = c("dose"),
    validation_fold = 5
  ) |>
  predict_survival(
    landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
    horizons = seq(from = 2 * 365.25, to = 6 * 365.25, by = 365.25),
    method = "coxph",
    type = "survival",
    dynamic_covariates = c("dose"),
    validation_fold = 5
  )
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> New names:
#> New names:
#> New names:
#> New names:
#> New names:
#> • `` -> `...10`
#> Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
#> Ran out of iterations and did not converge
#> New names:
#> New names:
#> New names:
#> New names:
#> New names:
#> • `` -> `...10`
```

We can also use [`summary()`](https://rdrr.io/r/base/summary.html) to
print parameter estimates. Note that this time the sample size is
smaller, because we have held out part of the original dataset for
validation.

``` r
summary(landmarking_object,
        type = "longitudinal",
        landmark = 365.25,
        dynamic_covariate = "dose")
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: value ~ treat + age + gender + learn.dis + (time | id)
#>    Data: dataframe
#> REML criterion at convergence: 1785.048
#> Random effects:
#>  Groups   Name        Std.Dev. Corr 
#>  id       (Intercept) 0.744917      
#>           time        0.003186 -0.16
#>  Residual             0.337911      
#> Number of obs: 862, groups:  id, 343
#> Fixed Effects:
#>  (Intercept)      treatLTG           age       genderM  learn.disYes  
#>    1.8730975    -0.1033299     0.0005706     0.1430885    -0.2732108  
#> optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 2 lme4 warnings
```

``` r
summary(landmarking_object,
        type = "survival",
        landmark = 365.25,
        horizon = 730.5)
#> Call:
#> survival::coxph(formula = formula, data = x@survival_datasets[[paste0(landmarks, 
#>     "-", horizons)]], model = TRUE, x = TRUE)
#> 
#>                   coef exp(coef)  se(coef)      z       p
#> treatLTG      0.050993  1.052316  0.214191  0.238 0.81182
#> age          -0.016284  0.983848  0.006482 -2.512 0.01199
#> genderM      -0.189162  0.827652  0.214458 -0.882 0.37775
#> learn.disYes -0.430140  0.650418  0.438529 -0.981 0.32666
#> dose          0.275454  1.317129  0.085492  3.222 0.00127
#> 
#> Likelihood ratio test=16.25  on 5 df, p=0.006154
#> n= 346, number of events= 90
```

Here are the in-sample performance metrics:

``` r
performance_metrics(
  landmarking_object,
  landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
  horizons = seq(from = 2 * 365.25, to = 6 * 365.25, by = 365.25),
  auc_t = TRUE
)
#>                landmark horizon    cindex     brier       AUCt
#> 365.25-730.5     365.25  730.50 0.6225997 0.5262495 0.93876976
#> 730.5-1095.75    730.50 1095.75 0.6432658 0.6500848 0.89580939
#> 1095.75-1461    1095.75 1461.00 0.6918492 0.6673227 0.84557518
#> 1461-1826.25    1461.00 1826.25 0.7616765 0.7284649 0.59756103
#> 1826.25-2191.5  1826.25 2191.50 0.9875992 0.9340199 0.07692308
```

Out-of-sample performance metrics can be obtained by specifying
`train = FALSE`:

``` r
performance_metrics(
  landmarking_object,
  landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
  horizons = seq(from = 2 * 365.25, to = 6 * 365.25, by = 365.25),
  auc_t = TRUE,
  train = FALSE
)
#> Warning in max(Y[status == 1]): no non-missing arguments to max; returning -Inf
#> Warning in max(Y[status == 1]): no non-missing arguments to max; returning -Inf
#>                landmark horizon    cindex     brier      AUCt
#> 365.25-730.5     365.25  730.50 0.5113732 0.5906314 0.9646766
#> 730.5-1095.75    730.50 1095.75 0.5475362 0.6748556 0.9823438
#> 1095.75-1461    1095.75 1461.00 0.1200000 0.7517810 1.0000000
#> 1461-1826.25    1461.00 1826.25       NaN 0.8027774        NA
#> 1826.25-2191.5  1826.25 2191.50       NaN 0.9843228        NA
```

## Cross-validation

Now, we can embed the above pipeline in a for loop to perform
cross-validation:

``` r
landmarking_object <- LandmarkAnalysis(
  data_static = static,
  data_dynamic = dynamic,
  event_indicator = "with.status",
  ids = "id",
  event_time = "with.time",
  times = "time",
  measurements = "value",
  K = 5
)

landmarking_object <- landmarking_object |>
  compute_risk_sets(
    landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25)
  )
```

``` r
metrics <- list()
for (k in 1:5) {
  metrics[[k]] <- landmarking_object |>
    fit_longitudinal(
      landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
      method = "lme4",
      formula = value ~ treat + age + gender + learn.dis + (time | id),
      dynamic_covariates = c("dose"),
      validation_fold = k
    ) |>
    predict_longitudinal(
      landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
      method = "lme4",
      allow.new.levels = TRUE,
      dynamic_covariates = c("dose"),
      validation_fold = k
    ) |>
    fit_survival(
      formula = Surv(event_time, event_status) ~
        treat + age + gender + learn.dis + dose,
      landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
      horizons = seq(from = 2 * 365.25, to = 6 * 365.25, by = 365.25),
      method = "coxph",
      dynamic_covariates = c("dose"),
      validation_fold = k
    ) |>
    predict_survival(
      landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
      horizons = seq(from = 2 * 365.25, to = 6 * 365.25, by = 365.25),
      method = "coxph",
      type = "survival",
      dynamic_covariates = c("dose"),
      validation_fold = k
    ) |>
    performance_metrics(
      landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
      horizons = seq(from = 2 * 365.25, to = 6 * 365.25, by = 365.25),
      auc_t = TRUE
    )
}
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> New names:
#> New names:
#> New names:
#> New names:
#> New names:
#> • `` -> `...10`
#> Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
#> Ran out of iterations and did not converge
#> New names:
#> New names:
#> New names:
#> New names:
#> New names:
#> • `` -> `...10`
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> New names:
#> New names:
#> New names:
#> New names:
#> New names:
#> • `` -> `...10`
#> Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
#> Ran out of iterations and did not converge
#> New names:
#> New names:
#> New names:
#> New names:
#> New names:
#> • `` -> `...10`
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> New names:
#> New names:
#> New names:
#> New names:
#> • `` -> `...10`
#> Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
#> Loglik converged before variable 4 ; coefficient may be infinite.
#> New names:
#> • `` -> `...10`
#> Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
#> Ran out of iterations and did not converge
#> New names:
#> New names:
#> New names:
#> New names:
#> New names:
#> • `` -> `...10`
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> New names:
#> New names:
#> New names:
#> New names:
#> New names:
#> • `` -> `...10`
#> Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
#> Ran out of iterations and did not converge
#> New names:
#> New names:
#> New names:
#> New names:
#> New names:
#> • `` -> `...10`
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> Warning in
#> predict.merMod(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]],
#> : unused arguments ignored
#> New names:
#> New names:
#> New names:
#> New names:
#> New names:
#> • `` -> `...10`
#> Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
#> Ran out of iterations and did not converge
#> New names:
#> New names:
#> New names:
#> New names:
#> New names:
#> • `` -> `...10`

metrics
#> [[1]]
#>                landmark horizon    cindex     brier      AUCt
#> 365.25-730.5     365.25  730.50 0.6336712 0.5476312 0.9343432
#> 730.5-1095.75    730.50 1095.75 0.6681284 0.6626849 0.8324416
#> 1095.75-1461    1095.75 1461.00 0.6973781 0.7012488 0.7499993
#> 1461-1826.25    1461.00 1826.25 0.7766896 0.7938033 0.4953194
#> 1826.25-2191.5  1826.25 2191.50 1.0000000 0.9846045 0.0000000
#> 
#> [[2]]
#>                landmark horizon    cindex     brier      AUCt
#> 365.25-730.5     365.25  730.50 0.6605406 0.5471592 0.9127225
#> 730.5-1095.75    730.50 1095.75 0.7141958 0.6657560 0.7919147
#> 1095.75-1461    1095.75 1461.00 0.7495549 0.7535546 0.7120984
#> 1461-1826.25    1461.00 1826.25 0.8210927 0.7853303 0.3812993
#> 1826.25-2191.5  1826.25 2191.50 1.0000000 0.9853555 0.0000000
#> 
#> [[3]]
#>                landmark horizon    cindex     brier      AUCt
#> 365.25-730.5     365.25  730.50 0.6212185 0.5485583 0.9260363
#> 730.5-1095.75    730.50 1095.75 0.6373537 0.6610433 0.9083269
#> 1095.75-1461    1095.75 1461.00 0.6475242 0.7164068 0.9113404
#> 1461-1826.25    1461.00 1826.25 0.7770478 0.7818805 0.6872417
#> 1826.25-2191.5  1826.25 2191.50 0.9875992 0.9300287 0.1153846
#> 
#> [[4]]
#>                landmark horizon    cindex     brier      AUCt
#> 365.25-730.5     365.25  730.50 0.5897559 0.5373165 0.9603579
#> 730.5-1095.75    730.50 1095.75 0.6387605 0.6681717 0.9376969
#> 1095.75-1461    1095.75 1461.00 0.7001645 0.7201016 0.7234967
#> 1461-1826.25    1461.00 1826.25 0.7562045 0.7458446 0.6171271
#> 1826.25-2191.5  1826.25 2191.50 0.9778120 0.9379297 0.0000000
#> 
#> [[5]]
#>                landmark horizon    cindex     brier       AUCt
#> 365.25-730.5     365.25  730.50 0.6059305 0.5527893 0.94672845
#> 730.5-1095.75    730.50 1095.75 0.6399187 0.6476502 0.85803018
#> 1095.75-1461    1095.75 1461.00 0.6417725 0.6504718 0.85764034
#> 1461-1826.25    1461.00 1826.25 0.7020792 0.7489997 0.67825788
#> 1826.25-2191.5  1826.25 2191.50 0.9870504 0.9262527 0.08333333
```
