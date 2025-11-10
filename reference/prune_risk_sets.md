# Prune a set of individuals from a risk set

Prune a set of individuals from a risk set

## Usage

``` r
prune_risk_sets(x, landmark, individuals)
```

## Arguments

- x:

  An object of class
  [`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md).

- landmark:

  a landmark time

- individuals:

  Vector of individuals to be pruned from

## Value

An object of class
[`LandmarkAnalysis`](https://vallejosgroup.github.io/landmaRk/reference/LandmarkAnalysis.md)
after having pruned the individuals indicated in `individuals` from the
risk set at landmark time `landmark`.
