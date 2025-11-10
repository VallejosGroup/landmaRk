# landmaRk

Time-to-event analysis using a wide array of longitudinal and survival
sub-models.

## Introduction

Time-to-event, or survival analysis, is used to analyse the time until
an *event of interest* occurs. Common events include hospitalisation,
equipment failure, or a prisoner reoffending. Whilst classic survival
methods assume model covariates are static, it is often the case that
longitudinal data related to the outcome of interest are collected. Two
main forms of survival analysis incorporating time-dependent covariates
exist, joint models and landmarking [¹](#fn1). This package focuses on
the latter.

For a set of landmark times, a survival model is fitted up to specified
horizon times. At landmark times, any time-dependent covariates must be
summarised. Most commonly, the last observation carried forward (LOCF)
approach is used. However, a more modern approach is to instead fit a
linear mixed effects model which accounts for observations being
measured with error [²](#fn2). However, any method which summarises
longitudinal observations can be used. Moreover, whilst landmarking
methods typically reply on Cox proportional hazards models, nearly any
survival model can also be used.

Whilst packages already exist which implement landmarking, these
packages implement specific longitudinal and survival models. The aim of
`landmaRk` is to support a wide array of longitudinal and survival
sub-models whilst providing a modular system allowing others to
incorporate their own models.

## Installation

We are planning to release the package on CRAN once the software is
mature. For now, you can get a CRAN-like experience by installing from
our r-universe

``` r
install.packages("landmaRk",
                 repos = c("https://vallejosgroup.r-universe.dev",
                           "https://cloud.r-project.org"))
```

Alternatively, the package can be built from source using `remotes`

``` r
# install.packages("remotes")
remotes::install_github("vallejosgroup/landmaRk", build_vignettes = TRUE)
```

## Getting started

We recommend starting with the `landmaRk` vignette, which provides an
overview of the package and how to use it. You can access the vignette
in R by calling

``` r
vignette("landmaRk")
```

Alternatively, you can view the vignette
[online](https://vallejosgroup.github.io/landmaRk/articles/landmaRk.html).

## References

------------------------------------------------------------------------

1.  Rizopoulos D, Molenberghs G, Lesaffre EMEH. Dynamic predictions with
    time-dependent covariates in survival analysis using joint modeling
    and landmarking. Biometrical Journal. 2017;59(6):1261-1276. doi:
    [10.1002/bimj.201600238](https://doi.org/10.1002/bimj.201600238)

2.  Paige E, Barrett J, Stevens D, et al. Landmark models for optimizing
    the use of repeated measurements of risk factors in electronic
    health records to predict future disease risk. American Journal of
    Epidemiology. 2018;187(7):1530-1538. doi:
    [10.1093/aje/kwy018](https://doi.org/10.1093/aje/kwy018)
