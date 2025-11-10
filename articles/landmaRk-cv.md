# Cross-validation with the landmaRk package

## Setup

In addition to the `landmaRk` package, we will also use `tidyverse`.

``` r
set.seed(123)
library(landmaRk)
library(tidyverse)
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.4     ✔ readr     2.1.5
#> ✔ forcats   1.0.1     ✔ stringr   1.6.0
#> ✔ ggplot2   4.0.0     ✔ tibble    3.3.0
#> ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
#> ✔ purrr     1.2.0     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
library(survival)
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
#>                   coef exp(coef)  se(coef)      z     p
#> treatLTG      0.042801  1.043730  0.289702  0.148 0.883
#> age          -0.012326  0.987750  0.008503 -1.450 0.147
#> genderM      -0.003815  0.996192  0.295850 -0.013 0.990
#> learn.disYes -0.860808  0.422820  0.740392 -1.163 0.245
#> dose          0.325584  1.384839  0.113134  2.878 0.004
#> 
#> Likelihood ratio test=10.34  on 5 df, p=0.06621
#> n= 346, number of events= 49
```

Here are the in-sample performance metrics:

``` r
performance_metrics(
  landmarking_object,
  landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
  horizons = seq(from = 2 * 365.25, to = 6 * 365.25, by = 365.25),
  auc_t = TRUE
)
#>                landmark horizon    cindex     brier     AUCt1     AUCt2
#> 365.25-730.5     365.25  730.50 0.2250417 0.6490448 0.9906947 0.9835355
#> 730.5-1095.75    730.50 1095.75 0.4002502 0.7862853 0.9629473 0.9331092
#> 1095.75-1461    1095.75 1461.00 0.4107143 0.7749309 0.9375000 0.7669928
#> 1461-1826.25    1461.00 1826.25 0.4311774 0.7548344 1.0000000 0.8785837
#> 1826.25-2191.5  1826.25 2191.50 0.9411765 0.9340199 0.1000000 0.1153846
#>                    AUCt3     AUCt4      AUCt5      AUCt6      AUCt7      AUCt8
#> 365.25-730.5   0.9343474 0.9262169 0.92984774 0.90451648 0.89613490 0.84661698
#> 730.5-1095.75  0.9099835 0.8668686 0.76147979 0.76349325 0.76209135 0.65998586
#> 1095.75-1461   0.7678328 0.7722389 0.68207760 0.66225974 0.59515953 0.59179396
#> 1461-1826.25   0.8986478 0.7956535 0.65172219 0.68151100 0.64682663 0.56038561
#> 1826.25-2191.5 0.1363636 0.1428571 0.04273504 0.04273504 0.05128205 0.05494505
#>                     AUCt9     AUCt10
#> 365.25-730.5   0.81209213 0.78216204
#> 730.5-1095.75  0.64683809 0.61057758
#> 1095.75-1461   0.57390978 0.57576113
#> 1461-1826.25   0.46876576 0.46026476
#> 1826.25-2191.5 0.06410256 0.06993007
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
#>                landmark horizon    cindex     brier     AUCt1     AUCt2
#> 365.25-730.5     365.25  730.50 0.1091234 0.7030735 1.0000000 1.0000000
#> 730.5-1095.75    730.50 1095.75 0.2583333 0.7902356 0.9565217 0.8701926
#> 1095.75-1461    1095.75 1461.00 0.0000000 0.8585820        NA        NA
#> 1461-1826.25    1461.00 1826.25       NaN 0.8296894        NA        NA
#> 1826.25-2191.5  1826.25 2191.50       NaN 0.9843228        NA        NA
#>                    AUCt3     AUCt4     AUCt5     AUCt6     AUCt7     AUCt8
#> 365.25-730.5   1.0000000 1.0000000 0.9855764 0.9925822 0.9969145 0.9969145
#> 730.5-1095.75  0.8229645 0.8321659 0.7002645 0.7074602 0.7174983 0.7174983
#> 1095.75-1461          NA 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
#> 1461-1826.25          NA        NA        NA        NA        NA        NA
#> 1826.25-2191.5        NA        NA        NA        NA        NA        NA
#>                    AUCt9    AUCt10
#> 365.25-730.5   0.9291742 0.9293744
#> 730.5-1095.75  0.7174983 0.7680128
#> 1095.75-1461   1.0000000 1.0000000
#> 1461-1826.25          NA        NA
#> 1826.25-2191.5        NA        NA
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
#> • `` -> `...10`
#> Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
#> Loglik converged before variable 4 ; coefficient may be infinite.
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
#> • `` -> `...10`
#> Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
#> Loglik converged before variable 4 ; coefficient may be infinite.
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
#>                landmark horizon    cindex     brier     AUCt1     AUCt2
#> 365.25-730.5     365.25  730.50 0.2576979 0.6750654 0.9847133 0.9692624
#> 730.5-1095.75    730.50 1095.75 0.4614756 0.7772294 0.8880067 0.7710943
#> 1095.75-1461    1095.75 1461.00 0.4652677 0.7796606 0.9414062 0.7370968
#> 1461-1826.25    1461.00 1826.25 0.5577342 0.8093310        NA 0.6836420
#> 1826.25-2191.5  1826.25 2191.50 1.0000000 0.9846045 0.0000000 0.0000000
#>                    AUCt3     AUCt4     AUCt5     AUCt6     AUCt7     AUCt8
#> 365.25-730.5   0.9307441 0.9212866 0.9096582 0.8849507 0.8518927 0.7956471
#> 730.5-1095.75  0.7058081 0.6779460 0.6092313 0.6114898 0.6172209 0.6000710
#> 1095.75-1461   0.7433333 0.7198520 0.6289192 0.6003492 0.5410070 0.5407021
#> 1461-1826.25   0.7061688 0.6248573 0.5220890 0.5650313 0.4550202 0.4482342
#> 1826.25-2191.5 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#>                    AUCt9    AUCt10
#> 365.25-730.5   0.7790505 0.7470220
#> 730.5-1095.75  0.5999108 0.5654821
#> 1095.75-1461   0.5184780 0.5251662
#> 1461-1826.25   0.3534164 0.3348095
#> 1826.25-2191.5 0.0000000 0.0000000
#> 
#> [[2]]
#>                landmark horizon    cindex     brier     AUCt1     AUCt2
#> 365.25-730.5     365.25  730.50 0.2798589 0.6660893 0.9881538 0.9735122
#> 730.5-1095.75    730.50 1095.75 0.4314570 0.7618801 0.9516195 0.8331193
#> 1095.75-1461    1095.75 1461.00 0.5960334 0.8469187 0.7203108 0.5583966
#> 1461-1826.25    1461.00 1826.25 0.5108514 0.7972720 0.8994652 0.8491005
#> 1826.25-2191.5  1826.25 2191.50 1.0000000 0.9853555        NA        NA
#>                    AUCt3     AUCt4     AUCt5     AUCt6     AUCt7     AUCt8
#> 365.25-730.5   0.8765133 0.8772083 0.8892168 0.8534006 0.8217506 0.7809371
#> 730.5-1095.75  0.7846446 0.7580944 0.6844112 0.6911200 0.6861711 0.6054836
#> 1095.75-1461   0.5624183 0.5618904 0.4258553 0.4276322 0.3973064 0.3952500
#> 1461-1826.25   0.7030205 0.7131008 0.5688396 0.5963637 0.5571901 0.4825371
#> 1826.25-2191.5        NA        NA 0.0000000 0.0000000 0.0000000 0.0000000
#>                    AUCt9    AUCt10
#> 365.25-730.5   0.7427626 0.7230432
#> 730.5-1095.75  0.5924479 0.5752228
#> 1095.75-1461   0.3907853 0.3895766
#> 1461-1826.25   0.3786107 0.3700142
#> 1826.25-2191.5 0.0000000 0.0000000
#> 
#> [[3]]
#>                landmark horizon    cindex     brier     AUCt1     AUCt2
#> 365.25-730.5     365.25  730.50 0.2713214 0.6716338 0.9889312 0.9878066
#> 730.5-1095.75    730.50 1095.75 0.4399887 0.7728541 0.8297670 0.7581583
#> 1095.75-1461    1095.75 1461.00 0.4491228 0.8065689 0.9365079 0.6472846
#> 1461-1826.25    1461.00 1826.25 0.3624679 0.8130292 0.9135802 0.8302186
#> 1826.25-2191.5  1826.25 2191.50 0.9215686 0.9300287 0.1290323 0.1538462
#>                    AUCt3     AUCt4     AUCt5     AUCt6      AUCt7     AUCt8
#> 365.25-730.5   0.8962888 0.8970480 0.9129100 0.8722543 0.84490634 0.7948612
#> 730.5-1095.75  0.7693931 0.7339311 0.6829594 0.6921019 0.69759775 0.6130315
#> 1095.75-1461   0.6553621 0.5860532 0.5249039 0.5844193 0.58829446 0.5412064
#> 1461-1826.25   0.8324390 0.8057394 0.8346940 0.8518402 0.79337907 0.6091133
#> 1826.25-2191.5 0.1818182 0.1904762 0.0678733 0.0678733 0.08241758 0.0887574
#>                    AUCt9    AUCt10
#> 365.25-730.5   0.7606236 0.7349016
#> 730.5-1095.75  0.6029425 0.5852143
#> 1095.75-1461   0.5488984 0.5566024
#> 1461-1826.25   0.6142147 0.6259838
#> 1826.25-2191.5 0.1048951 0.1153846
#> 
#> [[4]]
#>                landmark horizon    cindex     brier      AUCt1     AUCt2
#> 365.25-730.5     365.25  730.50 0.1439101 0.6483904 0.99571742 0.9931192
#> 730.5-1095.75    730.50 1095.75 0.3909853 0.7894662 0.94861100 0.8890465
#> 1095.75-1461    1095.75 1461.00 0.6751969 0.8353488 0.93650794 0.5426002
#> 1461-1826.25    1461.00 1826.25 0.4349593 0.7725270 1.00000000 0.8339512
#> 1826.25-2191.5  1826.25 2191.50 0.9600000 0.9379297 0.06666667 0.0800000
#>                    AUCt3      AUCt4     AUCt5     AUCt6     AUCt7     AUCt8
#> 365.25-730.5   0.9822635 0.98275289 0.9789553 0.9750663 0.9522820 0.9330249
#> 730.5-1095.75  0.8097079 0.78539769 0.6563541 0.6593042 0.6534115 0.6033175
#> 1095.75-1461   0.5394629 0.44873014 0.4509842 0.4520323 0.3720697 0.3515997
#> 1461-1826.25   0.8480321 0.77209800 0.6429344 0.6860274 0.6705238 0.5921040
#> 1826.25-2191.5 0.0500000 0.05263158 0.0000000 0.0000000 0.0000000 0.0000000
#>                    AUCt9    AUCt10
#> 365.25-730.5   0.9046571 0.8892080
#> 730.5-1095.75  0.5919417 0.6146831
#> 1095.75-1461   0.2930669 0.2968365
#> 1461-1826.25   0.5011351 0.5066887
#> 1826.25-2191.5 0.0000000 0.0000000
#> 
#> [[5]]
#>                landmark horizon    cindex     brier     AUCt1     AUCt2
#> 365.25-730.5     365.25  730.50 0.2279665 0.6823035 0.9897590 0.9857013
#> 730.5-1095.75    730.50 1095.75 0.4892277 0.7850246 0.9086824 0.7221125
#> 1095.75-1461    1095.75 1461.00 0.4051173 0.7451454 0.9791667 0.8062650
#> 1461-1826.25    1461.00 1826.25 0.4114833 0.7865477 0.9740260 0.8865393
#> 1826.25-2191.5  1826.25 2191.50 0.9361702 0.9262527 0.1071429 0.1304348
#>                    AUCt3     AUCt4     AUCt5     AUCt6      AUCt7      AUCt8
#> 365.25-730.5   0.9182847 0.9169552 0.9132184 0.8850731 0.88833970 0.85128205
#> 730.5-1095.75  0.6639545 0.6721947 0.6148379 0.6209510 0.63729112 0.55684347
#> 1095.75-1461   0.8100506 0.7655173 0.6754118 0.6693516 0.61024589 0.60889010
#> 1461-1826.25   0.8964118 0.8969084 0.7230070 0.7468498 0.74845052 0.57227809
#> 1826.25-2191.5 0.1578947 0.1666667 0.0468750 0.0468750 0.05357143 0.05357143
#>                     AUCt9    AUCt10
#> 365.25-730.5   0.82416651 0.7855345
#> 730.5-1095.75  0.55635200 0.5334360
#> 1095.75-1461   0.57881124 0.5820844
#> 1461-1826.25   0.50569049 0.4855763
#> 1826.25-2191.5 0.06818182 0.0750000
```
