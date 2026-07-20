initialise_visualisation_test_ <- function(K = 1) {
  epileptic_dfs <- split_wide_df(
    epileptic,
    ids = "id",
    times = "time",
    static = c(
      "with.time",
      "with.status",
      "treat",
      "age",
      "gender",
      "learn.dis"
    ),
    dynamic = c("dose"),
    measurement_name = "value"
  )

  LandmarkAnalysis(
    data_static = epileptic_dfs$df_static,
    data_dynamic = epileptic_dfs$df_dynamic,
    event_indicator = "with.status",
    ids = "id",
    event_time = "with.time",
    times = "time",
    measurements = "value",
    K = K
  ) |>
    compute_risk_sets(365.25)
}

# Colour-scale breaks double as the set of series shown in the legend
legend_breaks <- function(p) {
  p$scales$get_scales("colour")$breaks
}

test_that("plot() shows LOCF carried-forward value and survival curve", {
  withr::local_seed(123)
  x <- initialise_visualisation_test_() |>
    predict_longitudinal(
      landmarks = 365.25,
      method = "locf",
      dynamic_covariates = "dose"
    ) |>
    fit_survival(
      formula = Surv(event_time, event_status) ~ treat,
      landmarks = 365.25,
      horizons = 2 * 365.25,
      method = "coxph",
      dynamic_covariates = "dose"
    ) |>
    predict_survival(landmarks = 365.25, horizons = 2 * 365.25)

  one_id <- x@risk_sets[["365.25"]][1]
  p <- plot(x, id = one_id, landmark = 365.25, dynamic_covariate = "dose")

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, paste("Patient", one_id))
  expect_setequal(
    legend_breaks(p),
    c(
      "Observed measurements",
      "LOCF carried forward",
      "Predicted value (landmark)",
      "Survival probability"
    )
  )
  # No model-based trajectory curve for LOCF
  expect_false(any(vapply(
    p$layers,
    function(l) "average" %in% colnames(l$data),
    logical(1)
  )))
})

test_that("plot() shows population-average and individual trajectories for lme4", {
  withr::local_seed(123)
  x <- suppressWarnings(
    initialise_visualisation_test_() |>
      fit_longitudinal(
        landmarks = 365.25,
        method = "lme4",
        formula = value ~ treat + age + (time | id),
        dynamic_covariates = "dose"
      ) |>
      predict_longitudinal(
        landmarks = 365.25,
        method = "lme4",
        dynamic_covariates = "dose"
      ) |>
      fit_survival(
        formula = Surv(event_time, event_status) ~ treat,
        landmarks = 365.25,
        horizons = 2 * 365.25,
        method = "coxph",
        dynamic_covariates = "dose"
      ) |>
      predict_survival(landmarks = 365.25, horizons = 2 * 365.25)
  )

  one_id <- x@risk_sets[["365.25"]][1]
  p <- plot(x, id = one_id, landmark = 365.25, dynamic_covariate = "dose")

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, paste("Patient", one_id))
  expect_setequal(
    legend_breaks(p),
    c(
      "Observed measurements",
      "Population average trajectory",
      "Individual predicted trajectory",
      "Predicted value (landmark)",
      "Survival probability"
    )
  )
  # Every trajectory point should be at or before the landmark time
  trajectory_data <- p$layers[[
    which(vapply(
      p$layers,
      function(l) "average" %in% colnames(l$data),
      logical(1)
    ))[1]
  ]]$data
  expect_true(all(trajectory_data$time <= 365.25))
})

test_that("plot() shows cluster-average and individual trajectories for lcmm", {
  withr::local_seed(123)
  x <- suppressWarnings(
    initialise_visualisation_test_() |>
      fit_longitudinal(
        landmarks = 365.25,
        method = "lcmm",
        formula = value ~ treat + age + time,
        mixture = ~ treat + age + time,
        random = ~time,
        subject = "id",
        ng = 2,
        dynamic_covariates = "dose"
      ) |>
      predict_longitudinal(
        landmarks = 365.25,
        method = "lcmm",
        avg = FALSE,
        var.time = "time",
        dynamic_covariates = "dose"
      ) |>
      fit_survival(
        formula = Surv(event_time, event_status) ~ treat,
        landmarks = 365.25,
        horizons = 2 * 365.25,
        method = "coxph",
        dynamic_covariates = "dose"
      ) |>
      predict_survival(landmarks = 365.25, horizons = 2 * 365.25)
  )

  one_id <- x@risk_sets[["365.25"]][1]
  p <- suppressWarnings(
    plot(x, id = one_id, landmark = 365.25, dynamic_covariate = "dose")
  )

  expect_s3_class(p, "ggplot")
  expect_match(
    p$labels$title,
    paste0("^Patient ", one_id, " \\(most likely cluster: [12]\\)$")
  )
  expect_match(
    p$labels$subtitle,
    "^Cluster membership probabilities: cluster 1 = [01]\\.[0-9]{2}, cluster 2 = [01]\\.[0-9]{2}$"
  )
  breaks <- legend_breaks(p)
  # Both clusters' average trajectories are shown, not just the assigned one
  expect_true(all(
    c(
      "Cluster 1 average trajectory",
      "Cluster 2 average trajectory"
    ) %in%
      breaks
  ))
  expect_true("Individual predicted trajectory" %in% breaks)
  expect_true("Predicted value (landmark)" %in% breaks)
  expect_true("Survival probability" %in% breaks)

  # Check the trajectory data that plot() already computed (via its layers'
  # data), rather than calling the internal .lcmm_trajectory() helper a
  # second time in the same session: lcmm's prediction routines are backed
  # by Fortran code that is not always safe to invoke repeatedly, and doing
  # so has been observed to crash the R session
  avg_layer <- p$layers[[
    which(vapply(
      p$layers,
      function(l) "average" %in% colnames(l$data),
      logical(1)
    ))[1]
  ]]
  ind_layer <- p$layers[[
    which(vapply(
      p$layers,
      function(l) "individual" %in% colnames(l$data),
      logical(1)
    ))[1]
  ]]
  expect_setequal(levels(avg_layer$data$cluster), c("1", "2"))
  expect_equal(nrow(avg_layer$data), 200L)
  expect_equal(nrow(ind_layer$data), 100L)

  # Cross-check the probabilities reported in the subtitle against the
  # reported "most likely cluster" in the title
  probs <- as.numeric(regmatches(
    p$labels$subtitle,
    gregexpr("[01]\\.[0-9]{2}", p$labels$subtitle)
  )[[1]])
  cluster <- as.integer(regmatches(
    p$labels$title,
    regexpr("(?<=cluster: )[12]", p$labels$title, perl = TRUE)
  ))
  expect_equal(sum(probs), 1, tolerance = 1e-2)
  expect_equal(cluster, which.max(probs))
})

test_that("plot() shows out-of-sample (train = FALSE) results for LOCF", {
  withr::local_seed(123)
  x <- initialise_visualisation_test_(K = 5) |>
    predict_longitudinal(
      landmarks = 365.25,
      method = "locf",
      dynamic_covariates = "dose",
      validation_fold = 1
    ) |>
    fit_survival(
      formula = Surv(event_time, event_status) ~ treat,
      landmarks = 365.25,
      horizons = 2 * 365.25,
      method = "coxph",
      dynamic_covariates = "dose",
      validation_fold = 1
    ) |>
    predict_survival(
      landmarks = 365.25,
      horizons = 2 * 365.25,
      validation_fold = 1
    )

  test_ids <- x@cv_folds |> dplyr::filter(fold == 1) |> dplyr::pull(id)
  one_id <- intersect(x@risk_sets[["365.25"]], test_ids)[1]

  p <- plot(
    x,
    id = one_id,
    landmark = 365.25,
    dynamic_covariate = "dose",
    train = FALSE
  )

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, paste("Patient", one_id))
  expect_setequal(
    legend_breaks(p),
    c(
      "Observed measurements",
      "LOCF carried forward",
      "Predicted value (landmark)",
      "Survival probability"
    )
  )

  # The highlighted landmark prediction must come from the *test*-fold
  # predictions, not the training-fold ones
  pred_layer <- p$layers[[
    which(vapply(
      p$layers,
      function(l) identical(colnames(l$data), c("t", "p")),
      logical(1)
    ))[1]
  ]]
  expect_equal(
    pred_layer$data$p,
    unname(x@longitudinal_predictions_test[["365.25"]][["dose"]][
      as.character(one_id)
    ])
  )
  expect_false(
    one_id %in% names(x@longitudinal_predictions[["365.25"]][["dose"]])
  )
})

test_that("plot() shows out-of-sample (train = FALSE) trajectories for lme4", {
  withr::local_seed(123)
  x <- suppressWarnings(
    initialise_visualisation_test_(K = 5) |>
      fit_longitudinal(
        landmarks = 365.25,
        method = "lme4",
        formula = value ~ treat + age + (time | id),
        dynamic_covariates = "dose",
        validation_fold = 1
      ) |>
      predict_longitudinal(
        landmarks = 365.25,
        method = "lme4",
        dynamic_covariates = "dose",
        validation_fold = 1,
        allow.new.levels = TRUE
      ) |>
      fit_survival(
        formula = Surv(event_time, event_status) ~ treat,
        landmarks = 365.25,
        horizons = 2 * 365.25,
        method = "coxph",
        dynamic_covariates = "dose",
        validation_fold = 1
      ) |>
      predict_survival(
        landmarks = 365.25,
        horizons = 2 * 365.25,
        validation_fold = 1
      )
  )

  test_ids <- x@cv_folds |> dplyr::filter(fold == 1) |> dplyr::pull(id)
  one_id <- intersect(x@risk_sets[["365.25"]], test_ids)[1]

  p <- plot(
    x,
    id = one_id,
    landmark = 365.25,
    dynamic_covariate = "dose",
    train = FALSE
  )

  expect_s3_class(p, "ggplot")
  expect_setequal(
    legend_breaks(p),
    c(
      "Observed measurements",
      "Population average trajectory",
      "Individual predicted trajectory",
      "Predicted value (landmark)",
      "Survival probability"
    )
  )
  # Every trajectory point should be at or before the landmark time, same as
  # the in-sample case (the trajectory is computed from the model fit, which
  # does not depend on train/test membership of the plotted individual)
  trajectory_data <- p$layers[[
    which(vapply(
      p$layers,
      function(l) "average" %in% colnames(l$data),
      logical(1)
    ))[1]
  ]]$data
  expect_true(all(trajectory_data$time <= 365.25))
})

test_that("plot() errors informatively when out-of-sample predictions are unavailable", {
  withr::local_seed(123)
  x <- initialise_visualisation_test_() |>
    predict_longitudinal(
      landmarks = 365.25,
      method = "locf",
      dynamic_covariates = "dose"
    ) |>
    fit_survival(
      formula = Surv(event_time, event_status) ~ treat,
      landmarks = 365.25,
      horizons = 2 * 365.25,
      method = "coxph",
      dynamic_covariates = "dose"
    ) |>
    predict_survival(landmarks = 365.25, horizons = 2 * 365.25)

  one_id <- x@risk_sets[["365.25"]][1]

  expect_error(
    plot(
      x,
      id = one_id,
      landmark = 365.25,
      dynamic_covariate = "dose",
      train = FALSE
    ),
    paste(
      "No out-of-sample survival predictions found. Call",
      "predict_survival\\(\\) with validation_fold > 0."
    )
  )
})
