library(tidyverse)
epileptic <- joineR::epileptic
epileptic <- epileptic |>
  select(id, time, with.time, with.status, dose, treat, age, gender, learn.dis)

usethis::use_data(epileptic, overwrite = TRUE)
