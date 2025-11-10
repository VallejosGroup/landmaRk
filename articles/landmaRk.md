# Introduction to the landmaRk package

## Overview

The landmaRk package provides a framework for landmarking analysis of
time-to-event data with time-varying covariates. It allows users to
perform survival analysis using longitudinal data, fitting models to the
time-varying covariates, and then using these predictions in survival
models.

Given a time-to-event outcome \\T_i\\, a landmark time \\s\\ and a time
horizon \\s + w\\, the goal of a landmarking analysis is to estimate \\
\pi_i(s+w \vert s) = P(T_i \> s+w \vert T_i \ge s, X_i(s)), \\ where
\\X_i(s)\\ denotes a vector of covariates which may include time-varying
covariates.

A landmarking analysis of time-to-event data has two components:

- First, model the longitudinal trajectories of dynamic covariates,
  \\X_i(t)\\. Then use, the fitted model to make a prediction for
  \\X_i(s)\\, \\\hat{X}\_i(s)\\, at the landmark time, \\s\\.

- Second, fit a survival model of the time-to-event outcome,
  conditioning on the predicted value for \\X_i(s)\\, and potentially on
  additional static and dynamic covariates.

The `landmaRk` package allows users to use the following method for the
first component:

- Last Observation Carried Forward (LOCF), using the last measurement
  for \\X_i\\ recorded prior to \\s\\ as our prediction,
  \\\hat{X}\_i(s)\\.

- Linear mixed-effect (LME) model, as implemented in the `lme4` package.

- Latent class mixed model (LCMM), as implemented in the `lcmm` package.

For the second component, at present the `landmaRk` package supports Cox
proportional hazard models as implemented in the `survival`package.

Additionally, the `landmaRk` package provides a modular system allowing
making it possible to incorporate additional models both for the
longitudinal and the survival components

![Diagram of the landmaRk package
pipeline](../../../_temp/Library/landmaRk/diagram.svg)

Diagram of the landmaRk package pipeline

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
```

## Example: epileptic data

In this vignette, we use the dataset epileptic to perform landmarking
analysis of time-to-event data with time-varying covariates. Here is the
structure of the dataset.

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

The dataset contains the following variables:

- id: a unique patient identifier

- time: time when time-varying covariate “dose” was recorded

- with.time: time when the first of event or censoring happened

- with.status: indicates whether event (1) or censoring (0) occurred

- dose: a time-varying covariate

- treat, age, gender, learn.dis: static (baseline) covariates

## Initialising the landmarking analysis

First, we split the dataset into two, one containing static covariates,
event time and indicator of event/censoring, and another one containing
dynamic covariates. To that end, we use the function split_wide_df. That
function returns a named list with the following elements:

- Under the name df_static, a dataframe containing static covariates,
  event time and indicator of event/censoring.

- Under the name df_dynamic, a named list of dataframes, mapping dynamic
  covariates to dataframes in long format containing longitudinal
  measurement of the relevant dynamic covariate.

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

We can now create an object of class `LandmarkAnalysis`, using the
helper function of the same name.

``` r
landmarking_object <- LandmarkAnalysis(
  data_static = static,
  data_dynamic = dynamic,
  event_indicator = "with.status",
  ids = "id",
  event_time = "with.time",
  times = "time",
  measurements = "value"
)
```

Arguments to the helper function are the following:

- data_static and data_dynamic: the two datasets that were just created.

- event_indicator: name of the column that indicates the censoring
  indicator in the static dataset.

- dynamic_covariates: array column names in the dynamic dataset
  indicating time-varying covariates.

- ids: name of the column that identifies patients in both datasets.

- event_time: name of the column that identifies time of event/censoring
  in the static dataset.

## Baseline survival analysis (without time-varying covariates)

First, we perform a survival without time-varying covariates. We can use
this as a baseline to evaluate the performance of a subsequent landmark
analysis with such covariates. First step is to establish the landmark
times, and to work out the risk sets at each of those landmark times.

``` r
landmarking_object <- landmarking_object |>
  compute_risk_sets(landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25))

landmarking_object
#> Summary of LandmarkAnalysis Object:
#>   Landmarks: 365.25 730.5 1095.75 1461 1826.25 
#>   Number of observations: 605 
#>   Event indicator: with.status 
#>   Event time: with.time 
#>   Risk sets: 
#>     Landmark 365.25: 430 subjects
#>     Landmark 730.5: 270 subjects
#>     Landmark 1095.75: 168 subjects
#>     Landmark 1461: 111 subjects
#>     Landmark 1826.25: 47 subjects
```

Now we use the function `fit_survival` to fit the survival model. We
specify the following arguments:

- `landmarks`: Vector of landmark times at which the model will be
  fitted.

- `formula`: Two-sided formula object specifying the survival process.
  Because we are fitting the baseline model, we include static
  covariates only.

- `horizons`: Vector of time horizons up to when the model will be
  fitted. The vectors of landmarks and time horizons must be of the same
  length.

- `method`: Method for the survival component of the landmarking
  analysis. In this case, `"coxph"`.

- `dynamic_covariates`: Vector of names of the dynamic covariates to be
  used. In this baseline analysis, the vector is empty.

Then, we can make predictions using the function `predict_survival`.
Arguments `method`, `landmarks` and `horizons`are as above.

``` r
landmarking_object <- landmarking_object |>
  fit_survival(
    landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
    formula = Surv(event_time, event_status) ~ treat + age + gender + learn.dis,
    horizons = seq(from = 2 * 365.25, to = 6 * 365.25, by = 365.25),
    method = "coxph",
    dynamic_covariates = c()
  ) |>
  predict_survival(
    landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
    horizons = seq(from = 2 * 365.25, to = 6 * 365.25, by = 365.25),
    method = "coxph",
    type = "survival"
  )
#> Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
#> Ran out of iterations and did not converge
```

To display the results, one can use the method `summary`, specifying
`type = "survival"`, in addition to a landmark and a horizon.

``` r
summary(landmarking_object, type = "survival", landmark = 365.25, horizon = 730.5)
#> Call:
#> survival::coxph(formula = formula, data = x@survival_datasets[[paste0(landmarks, 
#>     "-", horizons)]], model = TRUE, x = TRUE)
#> 
#>                   coef exp(coef)  se(coef)      z     p
#> treatLTG      0.149026  1.160703  0.268802  0.554 0.579
#> age          -0.009630  0.990416  0.007368 -1.307 0.191
#> genderM       0.293059  1.340522  0.272513  1.075 0.282
#> learn.disYes -0.842526  0.430621  0.736509 -1.144 0.253
#> 
#> Likelihood ratio test=4.24  on 4 df, p=0.3746
#> n= 430, number of events= 57
```

Now the `performance_metrics` function can be used to calculate (for
now, in-sample) performance metrics.

``` r
performance_metrics(
  landmarking_object,
  landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
  horizons = seq(from = 2 * 365.25, to = 6 * 365.25, by = 365.25)
)
#>                landmark horizon    cindex     brier
#> 365.25-730.5     365.25  730.50 0.1298831 0.6611129
#> 730.5-1095.75    730.50 1095.75 0.3015631 0.7660603
#> 1095.75-1461    1095.75 1461.00 0.2382749 0.7905402
#> 1461-1826.25    1461.00 1826.25 0.3918367 0.7898793
#> 1826.25-2191.5  1826.25 2191.50 0.9384615 0.9321694
```

## Landmarking analysis with lme4 + coxph

Now we use the package lme4 to fit a linear mixed model of the
time-varying covariate, dose. This first step is followed by fitting a
Cox PH sub-model using the longitudinal predictions as covariates.

``` r
landmarking_object <- LandmarkAnalysis(
  data_static = static,
  data_dynamic = dynamic,
  event_indicator = "with.status",
  ids = "id",
  event_time = "with.time",
  times = "time",
  measurements = "value"
)

landmarking_object <- landmarking_object |>
  compute_risk_sets(
    landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25)
  )
```

Provided with the risk sets, now the pipeline has the following four
steps:

- Fit the longitudinal model with `fit_longitudinal`, specifying
  `method = "lme4"`, and the `formula` as it would be passed to
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html). One also
  needs to specify a vector of `landmarks` and a vector of dynamic
  covariates, `dynamic_covariates = c("dose")`.

- Make predictions with `predict_longitudinal`, specifying
  `method = "lme4"`.

- Fit the survival submodel with `fit_survival`, like in the baseline
  model section, but specifying the vector of dynamic covariates,
  `dynamic_covariates = c("dose")`.

- Make predictions with the survival submodel, like in the baseline
  model section.

``` r
landmarking_object <- landmarking_object |>
  fit_longitudinal(
    landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
    method = "lme4",
    formula = value ~ treat + age + gender + learn.dis + (time | id),
    dynamic_covariates = c("dose")
  ) |>
  predict_longitudinal(
    landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
    method = "lme4",
    allow.new.levels = TRUE,
    dynamic_covariates = c("dose")
  ) |>
  fit_survival(
    formula = Surv(event_time, event_status) ~
      treat + age + gender + learn.dis + dose,
    landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
    horizons = seq(from = 2 * 365.25, to = 6 * 365.25, by = 365.25),
    method = "coxph",
    dynamic_covariates = c("dose")
  ) |>
  predict_survival(
    landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
    horizons = seq(from = 2 * 365.25, to = 6 * 365.25, by = 365.25),
    method = "coxph",
    type = "survival"
  )
#> New names:
#> New names:
#> New names:
#> New names:
#> New names:
#> • `` -> `...10`
#> Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
#> Ran out of iterations and did not converge
```

As before, one can also use the function `summary` to display the
results.

``` r
summary(landmarking_object, type = "longitudinal", landmark = 365.25, dynamic_covariate = "dose")
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
```

``` r
summary(landmarking_object, type = "survival", landmark = 365.25, horizon = 730.5)
#> Call:
#> survival::coxph(formula = formula, data = x@survival_datasets[[paste0(landmarks, 
#>     "-", horizons)]], model = TRUE, x = TRUE)
#> 
#>                  coef exp(coef) se(coef)      z       p
#> treatLTG      0.12420   1.13224  0.26962  0.461 0.64506
#> age          -0.01015   0.98990  0.00776 -1.308 0.19098
#> genderM       0.20803   1.23126  0.27567  0.755 0.45045
#> learn.disYes -0.87223   0.41802  0.73770 -1.182 0.23706
#> dose          0.31780   1.37410  0.10911  2.913 0.00358
#> 
#> Likelihood ratio test=11.92  on 5 df, p=0.03587
#> n= 430, number of events= 57
```

Here are the performance metrics:

``` r
performance_metrics(
  landmarking_object,
  landmarks = seq(from = 365.25, to = 5 * 365.25, by = 365.25),
  horizons = seq(from = 2 * 365.25, to = 6 * 365.25, by = 365.25)
)
#>                landmark horizon    cindex     brier
#> 365.25-730.5     365.25  730.50 0.2241851 0.6672443
#> 730.5-1095.75    730.50 1095.75 0.4266108 0.7758282
#> 1095.75-1461    1095.75 1461.00 0.4021563 0.7967850
#> 1461-1826.25    1461.00 1826.25 0.3959184 0.7907983
#> 1826.25-2191.5  1826.25 2191.50 0.9384615 0.9445921
```

## Landmarking analysis with lcmm + coxph

Now we use the [lcmm
package](https://cecileproust-lima.github.io/lcmm/index.html) to fit a
latent class mixed model of the time-varying covariate, dose. This first
step is followed by fitting a Cox PH sub-model using the longitudinal
predictions as covariates.

The pipeline is identical to that of the lme4+coxph model. However, this
time one has to specify `method = "lcmm"` when calling
`fit_longitudinal`, in addition to the arguments that will be passed on
to
[`lcmm::hlme()`](https://cecileproust-lima.github.io/lcmm/reference/hlme.html),
namely

- `formula`: a two-sided formula of the fixed effects.

- `mixture`: a one-sided formula of the class-specific fixed effects.

- `random`: a one-sided formula specifying the random effects.

- `subject`: name of the column specifying the subject IDs.

- `ng`: the number of clusters in the LCMM model.

``` r
landmarking_object <- LandmarkAnalysis(
  data_static = static,
  data_dynamic = dynamic,
  event_indicator = "with.status",
  ids = "id",
  event_time = "with.time",
  times = "time",
  measurements = "value"
)
```

``` r
landmarking_object <- landmarking_object |>
  compute_risk_sets(seq(from = 365.25, to = 4 * 365.25, by = 365.25)) |>
  fit_longitudinal(
    landmarks = seq(from = 365.25, to = 4 * 365.25, by = 365.25),
    method = "lcmm",
    formula = value ~ treat + age + gender + learn.dis + time,
    mixture = ~ treat + age + gender + learn.dis + time,
    random = ~time,
    subject = "id",
    ng = 2,
    dynamic_covariates = c("dose")
  ) |>
  predict_longitudinal(
    landmarks = seq(from = 365.25, to = 4 * 365.25, by = 365.25),
    method = "lcmm",
    subject = "id",
    avg = FALSE,
    include_clusters = TRUE,
    var.time = "time",
    dynamic_covariates = c("dose")
  ) |>
  fit_survival(
    formula = Surv(event_time, event_status) ~
      treat + age + gender + learn.dis + dose,
    landmarks = seq(from = 365.25, to = 4 * 365.25, by = 365.25),
    horizons = seq(from = 2 * 365.25, to = 5 * 365.25, by = 365.25),
    method = "coxph",
    dynamic_covariates = c("dose"),
    include_clusters = TRUE
  ) |>
  predict_survival(
    landmarks = seq(from = 365.25, to = 4 * 365.25, by = 365.25),
    horizons = seq(from = 2 * 365.25, to = 5 * 365.25, by = 365.25),
    method = "coxph",
    dynamic_covariates = c("dose"),
    include_clusters = TRUE,
    type = "survival",
  )
#> Warning in
#> method(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]], :
#> Individuals 28, 389, 473, have not been used in LCMM model fitting. Imputing
#> values for those individuals
#> Warning in
#> method(x@longitudinal_fits[[as.character(landmarks)]][[dynamic_covariate]], :
#> Individuals 28, 389, 473, have not been used in LCMM model fitting. Imputing
#> values for those individuals
#> Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
#> Loglik converged before variable 5,10 ; coefficient may be infinite.
#> Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
#> Loglik converged before variable 10 ; coefficient may be infinite.
#> Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
#> Loglik converged before variable 5,9,10 ; coefficient may be infinite.
```

``` r
summary(landmarking_object,
        type = "longitudinal",
        landmark = 365.25,
        dynamic_covariate = "dose")
#> Heterogenous linear mixed model 
#>      fitted by maximum likelihood method 
#>  
#> lcmm::hlme(fixed = value ~ treat + age + gender + learn.dis + 
#>     time, mixture = ~treat + age + gender + learn.dis + time, 
#>     random = ~time, subject = "id", classmb = ~1, ng = 2, returndata = TRUE)
#>  
#> Statistical Model: 
#>      Dataset: NULL 
#>      Number of subjects: 427 
#>      Number of observations: 1074 
#>      Number of latent classes: 2 
#>      Number of parameters: 17  
#>  
#> Iteration process: 
#>      Convergence criteria satisfied 
#>      Number of iterations:  18 
#>      Convergence criteria: parameters= 1.1e-09 
#>                          : likelihood= 7e-07 
#>                          : second derivatives= 3.5e-13 
#>  
#> Goodness-of-fit statistics: 
#>      maximum log-likelihood: -1012.76  
#>      AIC: 2059.51  
#>      BIC: 2128.48  
#>  
#>  
#> Maximum Likelihood Estimates: 
#>  
#> Fixed effects in the class-membership model:
#> (the class of reference is the last class) 
#> 
#>                      coef      Se    Wald p-value
#> intercept class1  1.93379 0.20138   9.603 0.00000
#> 
#> Fixed effects in the longitudinal model:
#> 
#>                         coef      Se    Wald p-value
#> intercept class1     1.79296 0.10677  16.793 0.00000
#> intercept class2     1.50591 0.36890   4.082 0.00004
#> treatLTG class1     -0.10047 0.07110  -1.413 0.15763
#> treatLTG class2     -0.79701 0.32750  -2.434 0.01495
#> age class1          -0.00108 0.00189  -0.571 0.56818
#> age class2           0.00923 0.00783   1.179 0.23832
#> genderM class1       0.06565 0.07207   0.911 0.36233
#> genderM class2       0.76241 0.24669   3.091 0.00200
#> learn.disYes class1 -0.32565 0.15996  -2.036 0.04177
#> learn.disYes class2  0.56195 0.51759   1.086 0.27761
#> time class1          0.00089 0.00017   5.395 0.00000
#> time class2          0.00822 0.00050  16.283 0.00000
#> 
#> 
#> Variance-covariance matrix of the random-effects:
#>           intercept time
#> intercept   0.32775     
#> time        0.00009    0
#> 
#>                              coef      Se
#> Residual standard error:  0.36771 0.01301
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
#>                                 coef exp(coef)  se(coef)      z      p
#> treatLTG                   -0.049098  0.952088  0.299207 -0.164 0.8697
#> cluster_dose2                     NA        NA  0.000000     NA     NA
#> age                        -0.016270  0.983862  0.008495 -1.915 0.0555
#> genderM                     0.131476  1.140511  0.303117  0.434 0.6645
#> learn.disYes               -1.318162  0.267627  1.027134 -1.283 0.1994
#> dose                        0.411263  1.508723  0.257818  1.595 0.1107
#> treatLTG:cluster_dose2      0.670725  1.955655  0.920410  0.729 0.4662
#> cluster_dose2:age           0.042024  1.042919  0.020860  2.015 0.0439
#> cluster_dose2:genderM       0.497034  1.643838  0.848175  0.586 0.5579
#> cluster_dose2:learn.disYes  1.907475  6.736059  1.555057  1.227 0.2200
#> cluster_dose2:dose         -0.633257  0.530860  0.515536 -1.228 0.2193
#> 
#> Likelihood ratio test=12.71  on 10 df, p=0.2402
#> n= 430, number of events= 57
```

``` r
performance_metrics(
  landmarking_object,
  landmarks = seq(from = 365.25, to = 4 * 365.25, by = 365.25),
  horizons = seq(from = 2 * 365.25, to = 5 * 365.25, by = 365.25)
)
#>               landmark horizon    cindex     brier
#> 365.25-730.5    365.25  730.50 0.2652899 0.6752188
#> 730.5-1095.75   730.50 1095.75 0.5411742 0.7847448
#> 1095.75-1461   1095.75 1461.00 0.5283019 0.8091608
#> 1461-1826.25   1461.00 1826.25 0.6680272 0.8105012
```
