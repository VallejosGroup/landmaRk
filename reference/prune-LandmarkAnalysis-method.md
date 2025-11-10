# Prunes a landmark time from a [`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md), removing the risk set, longitudinal submodel and survival submodel from the object.

Prunes a landmark time from a
[`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md),
removing the risk set, longitudinal submodel and survival submodel from
the object.

## Usage

``` r
# S4 method for class 'LandmarkAnalysis'
prune(x, landmark = NULL)
```

## Arguments

- x:

  An object of class
  [`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md).

- landmark:

  A numeric indicating the landmark time.

## Value

An object of class
[`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md).
