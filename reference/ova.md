# Advanced ovarian cancer data

Data from two clinical trials comparing combination chemotherapies for
advanced ovarian cancer (Neijt et al. 1984, 1987). The data were
originally distributed in the archived dynpred package (as `ova`).

## Usage

``` r
data(ova)
```

## Format

A data frame with 358 rows, one per patient:

- `tyears`:

  numeric: time in years until death or last follow-up.

- `d`:

  survival status; `1 = `death, `0 = `censored.

- `Karn`:

  numeric: Karnofsky score, recoded from 6 (\\\le 60\\) to 10 (100).

- `Broders`:

  factor: Broders' histological grade, with levels `"1"`, `"2"`, `"3"`,
  `"4"`, `"unknown"`.

- `FIGO`:

  factor: FIGO stage, with levels `"III"` and `"IV"`.

- `Ascites`:

  factor: presence of ascites, with levels `"absent"`, `"present"`,
  `"unknown"`.

- `Diam`:

  factor: diameter of the residual tumour, with levels `"micr."`,
  `"<1cm"`, `"1-2cm"`, `"2-5cm"`, `">5cm"`.

- `id`:

  integer: patient identifier.

## Source

The archived dynpred package (GPL (\>= 2)), companion to van Houwelingen
& Putter (2012).

## References

van Houwelingen HC, Putter H (2012). *Dynamic Prediction in Clinical
Survival Analysis*. Chapman & Hall.

Neijt JP, ten Bokkel Huinink WW, van der Burg MEL, et al. Randomised
trial comparing two combination chemotherapy regimens (Hexa-CAF vs
CHAP-5) in advanced ovarian carcinoma. *Lancet*. 1984; **2**: 594-600.

Neijt JP, ten Bokkel Huinink WW, van der Burg MEL, et al. Randomized
trial comparing two combination chemotherapy regimens (CHAP-5 vs CP) in
advanced ovarian carcinoma. *J Clin Oncol*. 1987; **5**: 1157-1168.
