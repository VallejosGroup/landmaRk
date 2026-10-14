# Chronic myeloid leukaemia data: longitudinal white blood cell count

Repeated white blood cell count (WBC) measurements for the 190 patients
of [`cml`](https://vallejosgroup.github.io/landmaRk/reference/cml.md),
recorded every one to two months during the first stage of the Benelux
Chronic Myeloid Leukaemia (CML) trial (Kluin-Nelemans et al. 1998). The
data were originally distributed in the archived dynpred package (as
`wbc2`).

## Usage

``` r
data(cml_wbc)
```

## Format

A data frame with 5600 rows, in long format (one row per measurement):

- `id`:

  integer: patient identifier, shared with
  [`cml`](https://vallejosgroup.github.io/landmaRk/reference/cml.md).

- `time`:

  numeric: time in years from diagnosis at which the measurement was
  taken.

- `lwbc`:

  numeric: WBC, transformed and centered as
  \\{}^{10}\log(\mathrm{WBC}) - 0.95\\.

## Source

The archived dynpred package (GPL (\>= 2)), companion to van Houwelingen
& Putter (2012).

## References

van Houwelingen HC, Putter H (2012). *Dynamic Prediction in Clinical
Survival Analysis*. Chapman & Hall.

Kluin-Nelemans JC, Delannoy A, Louwagie A, le Cessie S, Hermans J, van
der Burgh JF, Hagemeijer AM, van den Berghe H, Benelux CML Study Group
(1998). Randomized study on hydroxyurea alone versus hydroxyurea
combined with low-dose interferon-alpha 2b for chronic myeloid leukemia.
*Blood*. 1998; **91**: 2713-2721.

## See also

[`cml`](https://vallejosgroup.github.io/landmaRk/reference/cml.md) for
the companion static covariates and survival outcome.
