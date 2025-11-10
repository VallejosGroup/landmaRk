# Concordance index for competing risks

Assess discriminative performance of predictions obtained from a
conventional or competing risks time-to-event model using time-dependent
concordance index.

## Usage

``` r
.CIndexCRisks(
  predictions,
  time,
  cens.code = 0,
  status,
  cause,
  tau,
  method = c("survival", "cifs")
)
```

## Arguments

- predictions:

  Numeric vector of model predictions.

- time:

  Numeric vector describing the time to the event of interest or
  censoring.

- cens.code:

  Value used to denote censoring in `status`. Defaults to 0.

- status:

  Vector of censoring status.

- cause:

  Event of interest.

- tau:

  Time c-index is evaluated.

- method:

  `'survival'` if the predictions are survival probabilities or `'cifs'`
  if they are cumulative incidence functions

## Value

Concordance index value.

## Details

Uses the proportion of correctly ordered risk pairs for the event \\k\\,
based on the predicted risk of the event up to time \\\tau\\.

\$\$C_k(\tau) = \frac{\sum\_{i=1}^N \sum\_{j=1}^N (A\_{ij} + B\_{ij})
\cdot Q\_{ij} \cdot N_i^k(\tau)}{\sum\_{i=1}^N \sum\_{j=1}^N (A\_{ij} +
B\_{ij}) \cdot N_i^k(\tau)}\$\$

A == risk ordering of patients, small time means patient 'i' at higher
risk than patient 'j' experiencing event of interest \\A\[i,j\] = 0\\
for tied event times.

B == risk ordering of patients, large time for patient 'i' means lower
risk than patient 'j' if not experienced the event of interest. Ties are
included in B

Q == the risk ordering of the subjects, i.e., is subject i assigned a
higher risk by the model than the subject j, for event \\E_k\\ until
time \\t\\. \\Q\[i,j\] = 0\\ for tied predictions.

N_t == number of subjects with survival time \< time point and
experience event of interest Tied event times are included

## References

Ahuja K, Schaar M van der. Joint Concordance Index. Published online
August 17, 2019.
[doi:10.48550/arXiv.1810.11207](https://doi.org/10.48550/arXiv.1810.11207)
