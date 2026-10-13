# Chronic Myeloid Leukaemia data. Distributed in the archived 'dynpred'
# package (GPL (>= 2)) as `wbc1` (one row per patient: survival outcome and
# static risk factors) and `wbc2` (long format: repeated WBC measurements).
library(dplyr)

tmp_static <- tempfile(fileext = ".RData")
download.file(
  "https://github.com/cran/dynpred/raw/master/data/wbc1.RData",
  tmp_static
)
load(tmp_static)

tmp_dynamic <- tempfile(fileext = ".RData")
download.file(
  "https://github.com/cran/dynpred/raw/master/data/wbc2.RData",
  tmp_dynamic
)
load(tmp_dynamic)

cml <- wbc1 |>
  rename(id = patnr) |>
  as.data.frame()

cml_wbc <- wbc2 |>
  rename(id = patnr, time = tyears) |>
  as.data.frame()
attributes(cml_wbc)[c("variable.labels", "codepage")] <- NULL

usethis::use_data(cml, overwrite = TRUE)
usethis::use_data(cml_wbc, overwrite = TRUE)
