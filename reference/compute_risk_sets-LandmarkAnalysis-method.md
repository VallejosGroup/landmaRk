# Compute the list of individuals at risk at landmark times

Compute the list of individuals at risk at landmark times

## Usage

``` r
# S4 method for class 'LandmarkAnalysis'
compute_risk_sets(x, landmarks, .warn_when_less_than = 0, ...)
```

## Arguments

- x:

  An object of class
  [`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md).

- landmarks:

  Numeric vector of landmark times

- .warn_when_less_than:

  Integer indicating that a warning will be raised when the number of
  observations prior to a landmark time is less than that integer for
  certain individuals.

- ...:

  Additional arguments (not used)

## Value

An object of class
[`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md),
including desired risk sets for the relevant landmark times.

## Details

A risk set describes all subjects still at risk (i.e., not experienced
the event of interest or censored) at a given time. In
[`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md),
risk sets define which subjects should be included in the longitudinal
and survival sub-models for each landmark time.

The risk sets are stored in the `risk_sets` slot of the
[`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md)
object, where each risk set is a list of indices corresponding to the
subjects at risk at the respective landmark time.
