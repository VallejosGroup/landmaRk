set.seed(123)
library(landmaRk)
library(lcmm)
library(tidyverse)
library(prodlim)
library(JMbayes2)

data("aids")
aids$patient <- as.numeric(aids$patient)
str(aids)


# DF with Static covariates
aids_dfs <- split_wide_df(
  aids,
  ids = "patient",
  times = "obstime",
  static = c("Time", "death", "drug", "gender", "prevOI"),
  dynamic = c("CD4"),
  measurement_name = "value"
)
static <- aids_dfs$df_static
head(static)

# DF with Dynamic covariates
dynamic <- aids_dfs$df_dynamic
head(dynamic[["CD4"]])

## General purpose function to plot trajectories and predictions

myfun <- function(landmarking_object, i) {

  # Data extraction
  longdf <- landmarking_object@data_dynamic$CD4
  survdf <- landmarking_object@survival_datasets
  coxph_1 <- landmarking_object@survival_fits[[1]]
  coxph_2 <- landmarking_object@survival_fits[[2]]
  survpred <- landmarking_object@survival_predictions

  # Manual predictions
  survpred_lp_1 <- predict(coxph_1, newdata = survdf[[1]][survdf[[1]]$patient == i,], type = "lp")
  survpred_lp_2 <- predict(coxph_2, newdata = survdf[[2]][survdf[[2]]$patient == i,], type = "lp")
  survpred_risk_1 <- predict(coxph_1, newdata = survdf[[1]][survdf[[1]]$patient == i,], type = "risk")
  survpred_risk_2 <- predict(coxph_2, newdata = survdf[[2]][survdf[[2]]$patient == i,], type = "risk")
  survpred_survival_1 <- predict(coxph_1, newdata = survdf[[1]][survdf[[1]]$patient == i,], type = "survival")
  survpred_survival_2 <- predict(coxph_2, newdata = survdf[[2]][survdf[[2]]$patient == i,], type = "survival")

  # Estimate survival curves
  survfit_1 <- survival::survfit(coxph_1,
                                 newdata = survdf[[1]][survdf[[1]]$patient == i,])
  survfit_2 <- survival::survfit(coxph_2,
                                 newdata = survdf[[2]][survdf[[2]]$patient == i,])

  # Plot
  par(mfrow = c(1,3))
  with(longdf[longdf$patient == i,], plot(obstime, value, type = "b"))
  abline(v = c(6, 8), col = c("red", "blue"))
  plot(survfit_1, main = "Landmark at 6 months")
  abline(v = 12, col = "red")
  abline(h = survpred[[1]][which(survdf[[1]]$patient == i)], col = "red", lty = 2)
  plot(survfit_2, main = "Landmark at 8 months")
  abline(v = 12, col = "blue")
  abline(h = survpred[[2]][which(survdf[[2]]$patient == i)], col = "blue", lty = 2)

  cat("Predictions: \n")
  cat(paste("Landmark 6:", survpred[[1]][which(survdf[[1]]$patient == i)], "\n"))
  cat(paste("Landmark 6 - lp:", survpred_lp_1, "\n"))
  cat(paste("Landmark 6 - risk:", survpred_risk_1, "\n"))
  cat(paste("Landmark 6 - survival:", survpred_survival_1, "\n"))
  cat(paste("Landmark 8:", survpred[[2]][which(survdf[[2]]$patient == i)], "\n"))
  cat(paste("Landmark 8 - lp:", survpred_lp_2, "\n"))
  cat(paste("Landmark 8 - risk:", survpred_risk_2, "\n"))
  cat(paste("Landmark 8 - survival:", survpred_survival_2, "\n"))
}


## No longitudinal covariates ####

landmarking_static <- LandmarkAnalysis(
  data_static = static,
  data_dynamic = dynamic,
  event_indicator = "death",
  ids = "patient",
  event_time = "Time",
  times = "obstime",
  measurements = "value"
)

landmarking_static <- landmarking_static |>
  compute_risk_sets(landmarks = c(6, 8))

landmarking_static <- landmarking_static |>
  fit_survival(
    formula = Surv(event_time, event_status) ~ drug,
    landmarks = c(6, 8),
    horizons = 12 + c(6, 8),
    method = "coxph",
    dynamic_covariates = c()
  ) |>
  predict_survival(
    landmarks = c(6, 8),
    horizons = 12 + c(6, 8),
    method = "coxph",
    type = "survival"
  )

### Plot trajectories and predictions for some individuals
myfun(landmarking_static, i = 454)

## LOCF ####

landmarking_locf <- LandmarkAnalysis(
  data_static = static,
  data_dynamic = dynamic,
  event_indicator = "death",
  ids = "patient",
  event_time = "Time",
  times = "obstime",
  measurements = "value"
)

landmarking_locf <- landmarking_locf |>
  compute_risk_sets(landmarks = c(6, 8))

landmarking_locf <- landmarking_locf |>
  predict_longitudinal(
    landmarks = c(6, 8),
    method = "locf",
    dynamic_covariates = c("CD4")
  ) |>
  fit_survival(
    formula = Surv(event_time, event_status) ~ drug,
    landmarks = c(6, 8),
    horizons = 12 + c(6, 8),
    method = "coxph",
    dynamic_covariates = c("CD4")
  ) |>
  predict_survival(
    landmarks = c(6, 8),
    horizons = 12 + c(6, 8),
    method = "coxph",
    dynamic_covariates = c("CD4"),
    type = "lp"
  )


myfun(landmarking_locf, i = 454)
