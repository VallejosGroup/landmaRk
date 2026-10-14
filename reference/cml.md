# Chronic myeloid leukaemia data: static covariates and survival outcome

Survival outcome and baseline (static) risk factors for 190 patients
from the Benelux Chronic Myeloid Leukaemia (CML) trial (Kluin-Nelemans
et al. 1998). The data were originally distributed in the archived
dynpred package (as `wbc1`).

## Usage

``` r
data(cml)
```

## Format

A data frame with 190 rows, one per patient:

- `id`:

  integer: patient identifier, shared with
  [`cml_wbc`](https://vallejosgroup.github.io/landmaRk/reference/cml_wbc.md).

- `tyears`:

  numeric: time in years until death or last follow-up.

- `d`:

  survival status; `1 = `death, `0 = `censored.

- `sokal`:

  numeric: Sokal score, a clinical risk index based on spleen size,
  percentage of circulating blasts, platelet count and age at diagnosis.

- `age`:

  numeric: age at diagnosis (years).

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

[`cml_wbc`](https://vallejosgroup.github.io/landmaRk/reference/cml_wbc.md)
for the companion longitudinal biomarker data.
