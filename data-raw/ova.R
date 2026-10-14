# Ovarian cancer data.
# Distributed in the archived 'dynpred' package (GPL (>= 2)).
tmp <- tempfile(fileext = ".RData")
download.file(
  "https://github.com/cran/dynpred/raw/master/data/ova.RData",
  tmp
)
load(tmp)
usethis::use_data(ova, overwrite = TRUE)
