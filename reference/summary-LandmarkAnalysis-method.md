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
