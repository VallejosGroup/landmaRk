# Builds a LandmarkAnalysis object with risk sets computed at `landmarks`,
# shared by the tests below.
.gridsearch_test_landmarks <- function(landmarks) {
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
    measurements = "value"
  ) |>
    compute_risk_sets(landmarks = landmarks)
}

test_that(".fit_lcmm sequential grid search (rep > 1) converges", {
  library(lcmm)
  withr::local_seed(1)
  data(data_hlme, package = "lcmm")

  fit <- .fit_lcmm(
    Y ~ Time * X1,
    data = data_hlme,
    mixture = ~Time,
    random = ~Time,
    subject = "ID",
    ng = 2,
    rep = 3,
    classmb = ~ X2 + X3,
    maxiter = 30
  )
  expect_s3_class(fit, "hlme")
  expect_equal(fit$conv, 1)
})

test_that(".fit_lcmm parallel grid search (rep > 1, cl > 1) converges", {
  library(lcmm)
  withr::local_seed(1)
  data(data_hlme, package = "lcmm")

  fit <- .fit_lcmm(
    Y ~ Time * X1,
    data = data_hlme,
    mixture = ~Time,
    random = ~Time,
    subject = "ID",
    ng = 2,
    rep = 3,
    classmb = ~ X2 + X3,
    maxiter = 30,
    cl = 2
  )
  expect_s3_class(fit, "hlme")
  expect_equal(fit$conv, 1)
  # model_fit$call$data is fixed up to the real data frame, not the
  # temporary global-environment binding used to work around
  # lcmm::gridsearch()'s cluster-export mechanism
  expect_true(is.data.frame(fit$call$data))
  expect_length(
    ls(envir = globalenv(), pattern = "^\\.landmaRk_gridsearch_data_"),
    0
  )
})

test_that(".fit_lcmm requires lcmm to be attached for rep > 1 grid search", {
  was_attached <- "package:lcmm" %in% search()
  if (was_attached) {
    suppressWarnings(detach("package:lcmm", unload = FALSE, force = TRUE))
  }
  withr::defer({
    if (was_attached) suppressPackageStartupMessages(library(lcmm))
  })

  data(data_hlme, package = "lcmm")
  expect_error(
    .fit_lcmm(
      Y ~ Time * X1,
      data = data_hlme,
      mixture = ~Time,
      random = ~Time,
      subject = "ID",
      ng = 2,
      rep = 2,
      classmb = ~ X2 + X3,
      maxiter = 30
    ),
    "requires the lcmm package to be attached"
  )
})

test_that("fit_longitudinal rejects cores > 1 combined with `cl`", {
  # lcmm::gridsearch()'s own (PSOCK) parallelism, nested inside landmaRk's
  # per-landmark (FORK) parallelism, causes sibling forked workers to race
  # for the same port and fail; fit_longitudinal() must reject this
  # combination rather than let it fail unpredictably.
  library(lcmm)
  x <- .gridsearch_test_landmarks(c(365.25, 730.5))

  expect_error(
    x |>
      fit_longitudinal(
        landmarks = c(365.25, 730.5),
        method = "lcmm",
        formula = value ~ treat + age + gender + learn.dis + time,
        mixture = ~ treat + age + gender + learn.dis + time,
        random = ~time,
        subject = "id",
        ng = 2,
        rep = 3,
        cores = 2,
        cl = 2,
        dynamic_covariates = "dose"
      ),
    "is not supported"
  )
})

# A lightweight custom fitting "method" (as supported by fit_longitudinal()'s
# `method` argument) that does no real model fitting -- it just records which
# OS process handled a given landmark. This isolates the `cores` (per-landmark)
# axis of parallelism from lcmm entirely, so the following two tests are fast
# and their PID evidence is not confounded by lcmm's own behaviour.
.pid_method <- function(formula, data, ...) {
  list(pid = Sys.getpid())
}

test_that("fit_longitudinal actually uses 2 distinct processes when cores = 2", {
  # cores > 1 relies on a FORK cluster, which does not exist on Windows;
  # fit_longitudinal() already falls back to single-threaded execution
  # there (see .supports_parallel()), so this assertion does not hold.
  skip_on_os("windows")

  landmarks <- seq(365.25, 4 * 365.25, by = 365.25)
  x <- .gridsearch_test_landmarks(landmarks) |>
    fit_longitudinal(
      landmarks = landmarks,
      method = .pid_method,
      formula = value ~ treat,
      dynamic_covariates = "dose",
      cores = 2
    )

  pids <- vapply(
    as.character(landmarks),
    function(lm) x@longitudinal_fits[[lm]][["dose"]]$pid,
    numeric(1)
  )
  cat("pids used across", length(landmarks), "landmarks (cores = 2):\n")
  print(pids)
  # 4 landmarks statically split across a 2-worker FORK cluster: exactly 2
  # distinct worker processes did the work, not just 2 processes requested
  expect_length(unique(pids), 2)
})

test_that("fit_longitudinal uses a single process across landmarks when cores = 1", {
  landmarks <- seq(365.25, 4 * 365.25, by = 365.25)
  x <- .gridsearch_test_landmarks(landmarks) |>
    fit_longitudinal(
      landmarks = landmarks,
      method = .pid_method,
      formula = value ~ treat,
      dynamic_covariates = "dose",
      cores = 1
    )

  pids <- vapply(
    as.character(landmarks),
    function(lm) x@longitudinal_fits[[lm]][["dose"]]$pid,
    numeric(1)
  )
  cat("pids used across", length(landmarks), "landmarks (cores = 1):\n")
  print(pids)
  expect_length(unique(pids), 1)
})

test_that("fit_longitudinal actually uses `cl`'s worker processes for lcmm grid search", {
  library(lcmm)
  withr::local_seed(123)

  # A cluster we create (and own) ourselves, rather than an integer, so we
  # can directly query its worker PIDs after fitting. lcmm::gridsearch()
  # never stops a caller-supplied cluster (only one it creates itself from
  # an integer `cl`), so it is still alive to query afterwards.
  my_cl <- parallel::makeCluster(2)
  withr::defer(parallel::stopCluster(my_cl))

  x <- .gridsearch_test_landmarks(365.25) |>
    fit_longitudinal(
      landmarks = 365.25,
      method = "lcmm",
      formula = value ~ treat + age + gender + learn.dis + time,
      mixture = ~ treat + age + gender + learn.dis + time,
      random = ~time,
      subject = "id",
      ng = 2,
      rep = 3,
      cl = my_cl,
      dynamic_covariates = "dose"
    )

  fit <- x@longitudinal_fits[["365.25"]][["dose"]]
  expect_s3_class(fit, "hlme")
  expect_equal(fit$conv, 1)

  worker_pids <- unlist(parallel::clusterCall(my_cl, Sys.getpid))
  cat("worker pids used by lcmm::gridsearch()'s `cl`:\n")
  print(worker_pids)
  # rep = 3 restarts statically split across our 2-worker cluster: both
  # workers actually computed at least one restart, not just 2 available
  expect_length(unique(worker_pids), 2)
})

test_that("fit_longitudinal's lcmm grid search runs single-process when `cl` is not set", {
  library(lcmm)
  withr::local_seed(123)

  x <- .gridsearch_test_landmarks(365.25) |>
    fit_longitudinal(
      landmarks = 365.25,
      method = "lcmm",
      formula = value ~ treat + age + gender + learn.dis + time,
      mixture = ~ treat + age + gender + learn.dis + time,
      random = ~time,
      subject = "id",
      ng = 2,
      rep = 3,
      dynamic_covariates = "dose"
    )

  fit <- x@longitudinal_fits[["365.25"]][["dose"]]
  expect_s3_class(fit, "hlme")
  expect_equal(fit$conv, 1)
  # No extra worker cluster was created for the grid search: the only
  # artifact a *parallel* grid search would leave behind (the temporary
  # global-environment data binding used to work around
  # lcmm::gridsearch()'s cluster-export mechanism) is absent, confirming
  # the plain, single-process lcmm::gridsearch() code path ran.
  leftover <- ls(envir = globalenv(), pattern = "^\\.landmaRk_gridsearch_data_")
  cat(
    "parallel-gridsearch artifacts left behind (expect none):",
    length(leftover),
    "\n"
  )
  expect_length(leftover, 0)
})
