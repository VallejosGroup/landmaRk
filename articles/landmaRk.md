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
library(lcmm)
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
library(prodlim)
```

## Example: `aids` data

In this vignette, we use the dataset `aids` to perform landmarking
analysis of time-to-event data with time-varying covariates. Here is the
structure of the dataset.

``` r
library(JMbayes2)
#> Loading required package: survival
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
#>     Landmark 6: 404 subjects
#>     Landmark 8: 379 subjects
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
    horizons = 12 + c(6, 8),
    method = "coxph",
    type = "lp"
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
#>           coef exp(coef) se(coef)    z      p
#> drugddI 0.3123    1.3665   0.1785 1.75 0.0802
#> 
#> Likelihood ratio test=3.08  on 1 df, p=0.07918
#> n= 404, number of events= 127
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
#>      landmark horizon   Brier(9) Brier(12) Brier(18)    AUC(3)    AUC(6)
#> 6-18        6      18 0.07900671 0.1651251 0.2342406 0.5437755 0.5539229
#> 8-20        8      20 0.10044734 0.1701394 0.2442496 0.5357143 0.5626498
#>        AUC(12)
#> 6-18 0.5022860
#> 8-20 0.4582525
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
    allow.new.levels = TRUE,
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
    horizons = 12 + c(6, 8),
    method = "coxph",
    type = "lp"
  )
#> New names:
#> New names:
#> • `` -> `...9`
```

As before, one can also use the function `summary` to display the
results.

``` r
summary(landmarking_object, type = "longitudinal", landmark = 6, dynamic_covariate = "CD4")
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: value ~ prevOI + obstime + (obstime | patient)
#>    Data: dataframe
#> REML criterion at convergence: 5308.102
#> Random effects:
#>  Groups   Name        Std.Dev. Corr
#>  patient  (Intercept) 3.9911       
#>           obstime     0.2094   0.00
#>  Residual             1.6903       
#> Number of obs: 1060, groups:  patient, 404
#> Fixed Effects:
#> (Intercept)   prevOIAIDS      obstime  
#>     10.4443      -4.5307      -0.1785  
#> optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings
```

``` r
summary(landmarking_object, type = "survival", landmark = 6, horizon = 18)
#> Call:
#> survival::coxph(formula = formula, data = x@survival_datasets[[paste0(landmarks, 
#>     "-", horizons)]], model = TRUE, x = TRUE)
#> 
#>           coef exp(coef) se(coef)    z      p
#> drugddI 0.3123    1.3665   0.1785 1.75 0.0802
#> 
#> Likelihood ratio test=3.08  on 1 df, p=0.07918
#> n= 404, number of events= 127
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
#>      landmark horizon   Brier(9) Brier(12) Brier(18)    AUC(3)    AUC(6)
#> 6-18        6      18 0.07900671 0.1651251 0.2342406 0.5437755 0.5539229
#> 8-20        8      20 0.10044734 0.1701394 0.2442496 0.5357143 0.5626498
#>        AUC(12)
#> 6-18 0.5022860
#> 8-20 0.4582525
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
    subject = "patient",
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
    method = "coxph",
    type = "lp",
    dynamic_covariates = c("CD4"),
    include_clusters = FALSE
  )
#> New names:
#> New names:
#> • `` -> `...9`
```

``` r
summary(landmarking_object, type = "longitudinal", landmark = 6, dynamic_covariate = "CD4")
#> Heterogenous linear mixed model 
#>      fitted by maximum likelihood method 
#>  
#> hlme(fixed = value ~ obstime + prevOI, mixture = ~obstime + prevOI, 
#>     random = ~obstime, subject = "patient", classmb = ~1, ng = 2, 
#>     nwg = TRUE, maxiter = maxiter, returndata = TRUE)
#>  
#> Statistical Model: 
#>      Dataset: NULL 
#>      Number of subjects: 404 
#>      Number of observations: 1060 
#>      Number of latent classes: 2 
#>      Number of parameters: 12  
#>  
#> Iteration process: 
#>      Convergence criteria satisfied 
#>      Number of iterations:  12 
#>      Convergence criteria: parameters= 7.1e-09 
#>                          : likelihood= 1.4e-08 
#>                          : second derivatives= 1e-12 
#>  
#> Goodness-of-fit statistics: 
#>      maximum log-likelihood: -2578.12  
#>      AIC: 5180.24  
#>      BIC: 5228.26  
#>  
#>  
#> Maximum Likelihood Estimates: 
#>  
#> Fixed effects in the class-membership model:
#> (the class of reference is the last class) 
#> 
#>                      coef      Se   Wald p-value
#> intercept class1  0.02129 0.17431  0.122 0.90279
#> 
#> Fixed effects in the longitudinal model:
#> 
#>                       coef      Se   Wald p-value
#> intercept class1   5.38287 0.31868 16.891 0.00000
#> intercept class2  13.51658 0.49286 27.425 0.00000
#> obstime class1    -0.16531 0.03393 -4.873 0.00000
#> obstime class2    -0.19030 0.04637 -4.104 0.00004
#> prevOIAIDS class1 -1.44437 0.31851 -4.535 0.00001
#> prevOIAIDS class2 -4.82749 0.68196 -7.079 0.00000
#> 
#> 
#> Variance-covariance matrix of the random-effects:
#>           intercept obstime
#> intercept  13.60174        
#> obstime    -0.24981 0.16909
#> 
#>                                     coef      Se
#> Proportional coefficient class1  0.33917 0.04167
#> Residual standard error:         1.54912 0.05467
```

``` r
summary(landmarking_object, type = "survival", landmark = 6, horizon = 18)
#> Call:
#> survival::coxph(formula = formula, data = x@survival_datasets[[paste0(landmarks, 
#>     "-", horizons)]], model = TRUE, x = TRUE)
#> 
#>           coef exp(coef) se(coef)    z      p
#> drugddI 0.3123    1.3665   0.1785 1.75 0.0802
#> 
#> Likelihood ratio test=3.08  on 1 df, p=0.07918
#> n= 404, number of events= 127
```

``` r
performance_metrics(
  landmarking_object,
  landmarks = c(6, 8),
  horizons = c(18, 20), 
  auc_t = TRUE, c_index = FALSE,
  h_times = c(3, 6, 12)
)
#>      landmark horizon   Brier(9) Brier(12) Brier(18)    AUC(3)    AUC(6)
#> 6-18        6      18 0.07900671 0.1651251 0.2342406 0.5437755 0.5539229
#> 8-20        8      20 0.10044734 0.1701394 0.2442496 0.5357143 0.5626498
#>        AUC(12)
#> 6-18 0.5022860
#> 8-20 0.4582525
```
