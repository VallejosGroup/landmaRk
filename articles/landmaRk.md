# Introduction to the landmaRk package

## Overview

The landmaRk package provides a framework for landmarking analysis of
time-to-event and longitudinal data. It allows users to perform dynamic
risk prediction for time-to-event outcomes whilst taking into account
longitudinal measurements (e.g. biomarkers measured over time).
Landmarking consists of a two-step framework. First, fitting models to
the time-varying covariates, and then using these predictions in
survival models.

Given a time-to-event outcome \\T_i\\, a landmark time \\s\\ and a time
horizon \\s + w\\, the goal of a landmarking analysis is to estimate \\
\pi_i(s+w \mid s) = P(T_i \> s+w \mid T_i \ge s, 𝒀_i(s)), \\ where
\\𝒀_i(s)\\ denotes a vector of covariates which may include time-varying
covariates.

A landmarking analysis of time-to-event data has two components:

- First, model the longitudinal trajectories of dynamic covariates,
  \\𝒀_i(t)\\. Then use, the fitted model to make a prediction for
  \\𝒀_i(s)\\, \\\hat{y}\_i(s)\\, at the landmark time, \\s\\.

- Second, fit a survival model of the time-to-event outcome,
  conditioning on the predicted value for \\𝒀_i(s)\\, and potentially on
  additional static and dynamic covariates.

The `landmaRk` package allows users to use the following method for the
first component:

- Last Observation Carried Forward (LOCF), using the last measurement
  for \\𝒀_i\\ recorded prior to \\s\\ as our prediction,
  \\\hat{y}\_i(s)\\.

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
#> Loading required package: survival
library(lcmm)
library(tidyverse)
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.2.1     ✔ readr     2.2.0
#> ✔ forcats   1.0.1     ✔ stringr   1.6.0
#> ✔ ggplot2   4.0.3     ✔ tibble    3.3.1
#> ✔ lubridate 1.9.5     ✔ tidyr     1.3.2
#> ✔ purrr     1.2.2
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
library(prodlim)
```

## Example: `aids` data

In this vignette, we use the dataset `aids` to perform landmarking
analysis of time-to-event data with time-varying covariates. Here is the
structure of the dataset.

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

The dataset contains the following variables:

- `patient`: a unique patient identifier

- `Time`: time when death or censoring was recorded

- `death`: indicates whether death (1) or censoring (0) occurred

- `obstime`: time when time-varying covariate `CD4` was recorded

- `CD4`: a time-varying covariate

- `drug`, `gender`, `prevOI`, `AZT`: static (baseline) covariates

## Initialising the landmarking analysis

First, we split the dataset into two, one containing static covariates,
event time and indicator of event/censoring, and another one containing
dynamic covariates. To that end, we use the function `split_wide_df`.
That function returns a named list with the following elements:

- Under the name `df_static`, a dataframe containing static covariates,
  event times and a binary indicator of event/censoring.

- Under the name `df_dynamic`, a named list of dataframes, mapping
  dynamic covariates to dataframes in long format containing
  longitudinal measurement of the relevant dynamic covariate.

The above split reduces data storage requirements, particularly for
large datasets with a large number of individuals or longitudinal
measurements. This is because static covariate values are stored only
once per individual, rather than repeatedly for each longitudinal
measurement.

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

We can now create an object of class `LandmarkAnalysis`, using the
helper function of the same name.

``` r

landmarking_object <- LandmarkAnalysis(
  data_static = static,
  data_dynamic = dynamic,
  event_indicator = "death",
  ids = "patient",
  event_time = "Time",
  times = "obstime",
  measurements = "value"
)
```

Arguments to the helper function are the following:

- `data_static` and `data_dynamic`: two datasets containing static and
  dynamic covariates, respectively (as created above using the
  `split_wide_df` function). Both datasets must contain a column with
  individual ids.

- `event_indicator`: name of the column that indicates the censoring
  indicator in `data_static`.

- `measurements`: name of the column in `data_dynamic` that contains the
  recorded values of the time-varying covariates (e.g., `"value"`).

- `times`: name of the column in `data_dynamic` that contains the
  measurement times associated with the time-varying covariates (e.g.,
  `"time"`).

- `ids`: name of the column that identifies patients in both datasets.

- `event_time`: name of the column that identifies time of
  event/censoring in the static dataset.

## Baseline survival analysis (without time-varying covariates)

First, we perform a survival without time-varying covariates. We can use
this as a baseline to evaluate the performance of a subsequent landmark
analysis with such covariates. First step is to establish the landmark
times, and to work out the risk sets at each of those landmark times.

``` r

landmarking_object <- landmarking_object |>
  compute_risk_sets(landmarks = c(6, 8))

landmarking_object
#> Summary of LandmarkAnalysis Object:
#>   Landmarks: 6 8 
#>   Number of observations: 467 
#>   Event indicator: death 
#>   Event time: Time 
#>   Risk sets: 
#>     Landmark 6: 403 subjects
#>     Landmark 8: 378 subjects
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
    formula = Surv(event_time, event_status) ~ drug,
    landmarks = c(6, 8),
    horizons = 12 + c(6, 8),
    method = "coxph",
    dynamic_covariates = c()
  ) |>
  predict_survival(
    landmarks = c(6, 8),
    horizons = 12 + c(6, 8)
  )
```

To display the results, one can use the method `summary`, specifying
`type = "survival"`, in addition to a landmark and a horizon.

``` r

summary(landmarking_object, type = "survival", landmark = 6, horizon = 18)
#> Call:
#> survival::coxph(formula = formula, data = x@survival_datasets[[paste0(landmarks, 
#>     "-", horizons)]], model = TRUE, x = TRUE)
#> 
#>           coef exp(coef) se(coef)     z      p
#> drugddI 0.3306    1.3918   0.1793 1.843 0.0653
#> 
#> Likelihood ratio test=3.42  on 1 df, p=0.06425
#> n= 403, number of events= 126
```

Now the `performance_metrics` function can be used to calculate (for
now, in-sample) performance metrics.

``` r

performance_metrics(
  landmarking_object,
  landmarks = c(6, 8),
  horizons = c(18, 20),
  auc_t = TRUE, c_index = FALSE,
  h_times = c(3, 6, 12)
)
#> Registered S3 method overwritten by 'cmprsk':
#>   method      from
#>   plot.cuminc lcmm
#>      landmark horizon   Brier(9) Brier(12) Brier(18)    AUC(9)   AUC(12)
#> 6-18        6      18 0.07706756 0.1638480 0.2338849 0.5521678 0.5573471
#> 8-20        8      20 0.09856622 0.1688013 0.2443032 0.5357143 0.5662291
#>        AUC(18)
#> 6-18 0.5040560
#> 8-20 0.4597691
```

## Landmarking analysis with lme4 + coxph

Now we use the package lme4 to fit a linear mixed model of the
time-varying covariate, CD4. This first step is followed by fitting a
Cox PH sub-model using the longitudinal predictions as covariates.

``` r

landmarking_object <- LandmarkAnalysis(
  data_static = static,
  data_dynamic = dynamic,
  event_indicator = "death",
  ids = "patient",
  event_time = "Time",
  times = "obstime",
  measurements = "value"
)
landmarking_object <- landmarking_object |>
  compute_risk_sets(landmarks = c(6, 8))
```

Provided with the risk sets, now the pipeline has the following four
steps:

- Fit the longitudinal model with `fit_longitudinal`, specifying
  `method = "lme4"`, and the `formula` as it would be passed to
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html). One also
  needs to specify a vector of `landmarks` and a vector of dynamic
  covariates, `dynamic_covariates = c("CD4")`.

- Make predictions with `predict_longitudinal`, specifying
  `method = "lme4"`.

- Fit the survival submodel with `fit_survival`, like in the baseline
  model section, but specifying the vector of dynamic covariates,
  `dynamic_covariates = c("CD4")`.

- Make predictions with the survival submodel, like in the baseline
  model section.

``` r

landmarking_object <- landmarking_object |>
  fit_longitudinal(
    landmarks = c(6, 8),
    method = "lme4",
    formula = value ~ prevOI + obstime + (obstime | patient),
    dynamic_covariates = c("CD4")
  ) |>
  predict_longitudinal(
    landmarks = c(6, 8),
    method = "lme4",
    dynamic_covariates = c("CD4")
  ) |>
  fit_survival(
    formula = Surv(event_time, event_status) ~ drug,
    landmarks = c(6, 8),
    horizons = 12 + c(6, 8),
    method = "coxph",
    dynamic_covariates = c("CD4")
  ) |>
  predict_survival(
    landmarks = c(6, 8),
    horizons = 12 + c(6, 8)
  )
```

As before, one can also use the function `summary` to display the
results.

``` r

summary(landmarking_object,
        type = "longitudinal",
        landmark = 6,
        dynamic_covariate = "CD4")
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: value ~ prevOI + obstime + (obstime | patient)
#>    Data: dataframe
#> REML criterion at convergence: 5298.972
#> Random effects:
#>  Groups   Name        Std.Dev. Corr 
#>  patient  (Intercept) 3.9915        
#>           obstime     0.2085   0.00 
#>  Residual             1.6926        
#> Number of obs: 1058, groups:  patient, 403
#> Fixed Effects:
#> (Intercept)   prevOIAIDS      obstime  
#>     10.4447      -4.5162      -0.1788
```

``` r

summary(landmarking_object, type = "survival", landmark = 6, horizon = 18)
#> Call:
#> survival::coxph(formula = formula, data = x@survival_datasets[[paste0(landmarks, 
#>     "-", horizons)]], model = TRUE, x = TRUE)
#> 
#>           coef exp(coef) se(coef)     z      p
#> drugddI 0.3306    1.3918   0.1793 1.843 0.0653
#> 
#> Likelihood ratio test=3.42  on 1 df, p=0.06425
#> n= 403, number of events= 126
```

Here are the performance metrics:

``` r

performance_metrics(
  landmarking_object,
  landmarks = c(6, 8),
  horizons = c(18, 20),
  auc_t = TRUE, c_index = FALSE,
  h_times = c(3, 6, 12)
)
#>      landmark horizon   Brier(9) Brier(12) Brier(18)    AUC(9)   AUC(12)
#> 6-18        6      18 0.07706756 0.1638480 0.2338849 0.5521678 0.5573471
#> 8-20        8      20 0.09856622 0.1688013 0.2443032 0.5357143 0.5662291
#>        AUC(18)
#> 6-18 0.5040560
#> 8-20 0.4597691
```

## Landmarking analysis with lcmm + coxph

Now we use the [lcmm
package](https://cecileproust-lima.github.io/lcmm/index.html) to fit a
latent class mixed model of the time-varying covariate, CD4. This first
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
  event_indicator = "death",
  ids = "patient",
  event_time = "Time",
  times = "obstime",
  measurements = "value"
)
```

``` r

landmarking_object <- landmarking_object |>
  compute_risk_sets(landmarks = c(6, 8)) |>
  fit_longitudinal(
    landmarks = c(6, 8),
    method = "lcmm",
    formula = value ~ obstime + prevOI,
    mixture = ~ obstime + prevOI,
    random = ~ obstime,
    subject = "patient",
    ng = 2,
    dynamic_covariates = c("CD4"),
    maxiter = 5000, rep = 25, nwg = TRUE
  ) |>
  predict_longitudinal(
    landmarks = c(6, 8),
    method = "lcmm",
    avg = TRUE,
    include_clusters = FALSE,
    var.time = "obstime",
    dynamic_covariates = c("CD4")
  ) |>
  fit_survival(
    formula = Surv(event_time, event_status) ~ drug,
    landmarks = c(6, 8),
    horizons = 12 + c(6, 8),
    method = "coxph",
    dynamic_covariates = c("CD4"),
    include_clusters = FALSE
  ) |>
  predict_survival(
    landmarks = c(6, 8),
    horizons = 12 + c(6, 8),
    dynamic_covariates = c("CD4"),
    include_clusters = FALSE
  )
```

``` r

summary(landmarking_object,
        type = "longitudinal",
        landmark = 6,
        dynamic_covariate = "CD4")
#> Heterogenous linear mixed model 
#>      fitted by maximum likelihood method 
#>  
#> hlme(fixed = value ~ obstime + prevOI, mixture = ~obstime + prevOI, 
#>     random = ~obstime, subject = "patient", classmb = ~1, ng = 2, 
#>     nwg = TRUE, maxiter = 24000, returndata = TRUE)
#>  
#> Statistical Model: 
#>      Dataset: NULL 
#>      Number of subjects: 403 
#>      Number of observations: 1058 
#>      Number of latent classes: 2 
#>      Number of parameters: 12  
#>  
#> Iteration process: 
#>      Convergence criteria satisfied 
#>      Number of iterations:  1 
#>      Convergence criteria: parameters= 1.1e-10 
#>                          : likelihood= 5e-10 
#>                          : second derivatives= 8.3e-11 
#>  
#> Goodness-of-fit statistics: 
#>      maximum log-likelihood: -2573.92  
#>      AIC: 5171.84  
#>      BIC: 5219.82  
#>  
#>  
#> Maximum Likelihood Estimates: 
#>  
#> Fixed effects in the class-membership model:
#> (the class of reference is the last class) 
#> 
#>                      coef      Se   Wald p-value
#> intercept class1  0.02034 0.17564  0.116 0.90782
#> 
#> Fixed effects in the longitudinal model:
#> 
#>                       coef      Se   Wald p-value
#> intercept class1   5.38687 0.31907 16.883 0.00000
#> intercept class2  13.51821 0.49332 27.402 0.00000
#> obstime class1    -0.16632 0.03397 -4.896 0.00000
#> obstime class2    -0.19019 0.04641 -4.098 0.00004
#> prevOIAIDS class1 -1.43232 0.31900 -4.490 0.00001
#> prevOIAIDS class2 -4.81769 0.68299 -7.054 0.00000
#> 
#> 
#> Variance-covariance matrix of the random-effects:
#>           intercept obstime
#> intercept  13.59457        
#> obstime    -0.24938 0.16877
#> 
#>                                     coef      Se
#> Proportional coefficient class1  0.33971 0.04181
#> Residual standard error:         1.55098 0.05479
```

``` r

summary(landmarking_object, type = "survival", landmark = 6, horizon = 18)
#> Call:
#> survival::coxph(formula = formula, data = x@survival_datasets[[paste0(landmarks, 
#>     "-", horizons)]], model = TRUE, x = TRUE)
#> 
#>           coef exp(coef) se(coef)     z      p
#> drugddI 0.3306    1.3918   0.1793 1.843 0.0653
#> 
#> Likelihood ratio test=3.42  on 1 df, p=0.06425
#> n= 403, number of events= 126
```

``` r

performance_metrics(
  landmarking_object,
  landmarks = c(6, 8),
  horizons = c(18, 20),
  auc_t = TRUE, c_index = FALSE,
  h_times = c(3, 6, 12)
)
#>      landmark horizon   Brier(9) Brier(12) Brier(18)    AUC(9)   AUC(12)
#> 6-18        6      18 0.07706756 0.1638480 0.2338849 0.5521678 0.5573471
#> 8-20        8      20 0.09856622 0.1688013 0.2443032 0.5357143 0.5662291
#>        AUC(18)
#> 6-18 0.5040560
#> 8-20 0.4597691
```
