# Build a fine time grid for one individual, holding their static covariate
# profile fixed and varying time from (at or before) 0 through to the
# landmark time. Used as `newdata` for model-based trajectory predictions.
.trajectory_newdata <- function(
  x,
  dynamic_covariate,
  id,
  landmark,
  n_pts = 100L
) {
  obs_times <- x@data_dynamic[[dynamic_covariate]][[x@times]]
  t_min <- min(0, suppressWarnings(min(obs_times, na.rm = TRUE)))
  grid <- seq(t_min, landmark, length.out = n_pts)

  # Base-R subsetting (rather than dplyr::filter()) sidesteps data-masking
  # ambiguity between the `id` argument and any `id`-named column
  static_row <- x@data_static[x@data_static[[x@ids]] == id, , drop = FALSE]
  if (nrow(static_row) != 1L) {
    stop(
      "Individual ",
      id,
      " must appear exactly once in @data_static (found ",
      nrow(static_row),
      ")."
    )
  }

  newdata <- static_row[rep(1L, length(grid)), , drop = FALSE]
  newdata[[x@times]] <- grid
  newdata[[x@ids]] <- id
  rownames(newdata) <- NULL
  newdata
}

# Empirical-Bayes (BLUP) random-effect contribution, for a single
# out-of-sample subject, evaluated over every row of `newdata_grid` (e.g. a
# fine time grid for trajectory plotting). Mirrors the mixed-model-equations
# calculation that .predict_lme4() (R/longitudinal_lme4.R) uses to produce
# out-of-sample landmark predictions, specialised to a single subject whose
# random effects are estimated from their own repeated measurements
# (`newdata_long`), so that the plotted trajectory agrees with the plotted
# landmark prediction for out-of-sample individuals.
.lme4_blup_trajectory <- function(fit, re_terms, newdata_grid, newdata_long) {
  sigma2 <- sigma(fit)^2
  beta_hat <- lme4::fixef(fit)
  X_long <- model.matrix(formula(fit, fixed.only = TRUE), newdata_long)
  attributes(X_long)$assign <- NULL
  attributes(X_long)$contrasts <- NULL

  Sigma_b <- as.matrix(unclass(lme4::VarCorr(fit)[[1]]))
  attributes(Sigma_b)$correlation <- NULL
  attributes(Sigma_b)$stddev <- NULL

  Z_long <- model.matrix(reformulate(deparse(re_terms[[2]])), newdata_long)
  response <- as.character(formula(fit, fixed.only = TRUE)[[2]])

  b_hat <- solve(crossprod(Z_long) / sigma2 + solve(Sigma_b)) %*%
    (t(Z_long) / sigma2) %*%
    (newdata_long[[response]] - as.vector(X_long %*% beta_hat))

  Z_grid <- model.matrix(reformulate(deparse(re_terms[[2]])), newdata_grid)
  as.vector(Z_grid %*% b_hat)
}

# Population-average (fixed-effects only) and individual-specific (BLUP)
# trajectory curves from a fitted lme4 model, for one individual, evaluated
# over a fine time grid up to the landmark time.
#
# When `id` was part of the data used to fit `fit`, lme4's own fitted random
# effects already reflect them, so a plain predict() call gives the correct
# subject-specific curve. When `id` is out-of-sample (not among the fitted
# grouping-factor levels), predict(..., allow.new.levels = TRUE) would
# silently zero out their random effect and return the population-average
# curve instead; .lme4_blup_trajectory() is used in that case to estimate
# their random effects from their own pre-landmark observations, matching
# the out-of-sample landmark prediction computed by .predict_lme4().
#
# Returns a list with:
#  - avg_curves: data.frame(time, cluster, average); a single "cluster"
#    (there is no latent-class structure for lme4)
#  - ind_curve: data.frame(time, individual)
#  - cluster / probs: NULL (no latent-class structure for lme4)
.lme4_trajectory <- function(
  fit,
  x,
  dynamic_covariate,
  id,
  landmark,
  n_pts = 100L
) {
  newdata <- .trajectory_newdata(x, dynamic_covariate, id, landmark, n_pts)
  marginal <- as.numeric(predict(
    fit,
    newdata = newdata,
    re.form = NA,
    allow.new.levels = TRUE
  ))

  bars <- lme4::findbars(formula(fit))
  is_fitted_id <- FALSE
  if (length(bars) == 1L) {
    grouping_var <- deparse(bars[[1]][[3]])
    fitted_ids <- rownames(lme4::ranef(fit)[[grouping_var]])
    is_fitted_id <- as.character(id) %in% fitted_ids
  }

  ind_pred <- if (is_fitted_id) {
    as.numeric(predict(fit, newdata = newdata, allow.new.levels = TRUE))
  } else if (length(bars) == 1L) {
    obs_data <- x@data_dynamic[[dynamic_covariate]]
    obs_data <- obs_data[
      obs_data[[x@ids]] == id & obs_data[[x@times]] <= landmark,
      ,
      drop = FALSE
    ]
    if (nrow(obs_data) == 0L) {
      marginal
    } else {
      newdata_long <- dplyr::left_join(obs_data, x@data_static, by = x@ids)
      marginal +
        .lme4_blup_trajectory(fit, bars[[1]], newdata, newdata_long)
    }
  } else {
    # No single random-effect grouping to estimate BLUPs from (either no
    # random effects, or more than one grouping); fall back to a plain
    # predict() call, as before
    as.numeric(predict(fit, newdata = newdata, allow.new.levels = TRUE))
  }

  list(
    avg_curves = data.frame(
      time = newdata[[x@times]],
      cluster = factor(1L),
      average = marginal
    ),
    ind_curve = data.frame(
      time = newdata[[x@times]],
      individual = ind_pred
    ),
    cluster = NULL,
    probs = NULL
  )
}

# Cluster-average (class-specific fixed-effects only) trajectory curves for
# *every* latent class, and an individual-specific trajectory (class-specific
# fixed effects plus predicted random effects, averaged across classes using
# the individual's posterior class-membership probabilities), from a fitted
# lcmm model, for one individual, evaluated over a fine time grid up to the
# landmark time.
#
# Returns a list with:
#  - avg_curves: data.frame(time, cluster, average), stacked over all ng
#    classes
#  - ind_curve: data.frame(time, individual): the posterior-probability-
#    weighted average of the class-specific subject-specific predictions
#  - cluster: the individual's most likely class (or class 1 for
#    single-class models / when no pre-landmark observations are available)
#  - probs: named numeric vector of posterior class-membership probabilities
#    (names are the class labels "1", "2", ...)
.lcmm_trajectory <- function(
  fit,
  x,
  dynamic_covariate,
  id,
  landmark,
  n_pts = 100L
) {
  hlme <- NULL # Global var (see .predict_lcmm())
  # .fit_lcmm() stores a namespace-qualified call (lcmm::hlme); lcmm's
  # predictY()/predictRE()/predictClass() reconstruct calls internally via
  # do.call() and error on that, so it must be unqualified first (mirrors
  # the same fix-up in .predict_lcmm())
  fit$call[[1]] <- rlang::expr(hlme)

  ng <- fit$ng
  newdata <- .trajectory_newdata(x, dynamic_covariate, id, landmark, n_pts)
  grid_times <- newdata[[x@times]]

  # Class-specific population-average predictions (one column per class)
  class_pred <- as.matrix(
    lcmm::predictY(fit, newdata = newdata, var.time = x@times)$pred
  )

  dyn_df <- x@data_dynamic[[dynamic_covariate]]
  obs_data <- dyn_df[
    dyn_df[[x@ids]] == id & dyn_df[[x@times]] <= landmark,
    ,
    drop = FALSE
  ]
  # predictClass()/predictRE() need every covariate in the model formulas
  # (e.g. mixture/classmb terms), not just the dynamic covariate itself
  if (nrow(obs_data) > 0L) {
    obs_data <- dplyr::left_join(obs_data, x@data_static, by = x@ids)
  }

  # Posterior class-membership probabilities for this individual
  if (ng == 1L) {
    probs <- stats::setNames(1, "1")
  } else if (nrow(obs_data) == 0L) {
    # No pre-landmark data to base a class assignment on: fall back to the
    # population-average class-membership probabilities (same imputation
    # idea as .predict_lcmm())
    prob_cols <- grep("^prob", colnames(fit$pprob), value = TRUE)
    probs <- stats::setNames(
      colMeans(fit$pprob[, prob_cols, drop = FALSE]),
      as.character(seq_len(ng))
    )
  } else {
    class_probs <- lcmm::predictClass(fit, newdata = obs_data, subject = x@ids)
    prob_cols <- grep("^prob", colnames(class_probs), value = TRUE)
    probs <- stats::setNames(
      as.numeric(class_probs[1L, prob_cols]),
      as.character(seq_len(ng))
    )
  }
  cluster <- as.integer(which.max(probs))

  # Individual (subject-specific) predictions require predicted random
  # effects, which in turn require the individual's own pre-landmark
  # observations; fall back to the (RE-free) class-average predictions
  # otherwise
  if (nrow(obs_data) > 0L) {
    predRE <- lcmm::predictRE(
      fit,
      obs_data,
      subject = x@ids,
      classpredRE = TRUE
    )
    ind_pred <- as.matrix(
      lcmm::predictY(
        fit,
        newdata = newdata,
        var.time = x@times,
        predRE = predRE
      )$pred
    )
  } else {
    ind_pred <- class_pred
  }
  # Average the class-specific individual predictions, weighted by the
  # individual's posterior class-membership probabilities
  individual <- as.numeric(ind_pred %*% probs)

  list(
    avg_curves = data.frame(
      time = rep(grid_times, ng),
      cluster = factor(rep(seq_len(ng), each = length(grid_times))),
      average = as.vector(class_pred)
    ),
    ind_curve = data.frame(
      time = grid_times,
      individual = individual
    ),
    cluster = cluster,
    probs = probs
  )
}

# Dispatch to the appropriate model-based trajectory helper, based on the
# class of the fitted longitudinal model. Returns NULL (with a warning) for
# unsupported model types or if trajectory computation fails, so that
# plot() can gracefully fall back to showing only the landmark-time
# prediction.
.model_trajectory <- function(
  fit,
  x,
  dynamic_covariate,
  id,
  landmark,
  n_pts = 100L
) {
  tryCatch(
    {
      if (is(fit, "hlme")) {
        .lcmm_trajectory(fit, x, dynamic_covariate, id, landmark, n_pts)
      } else if (is(fit, "merMod")) {
        .lme4_trajectory(fit, x, dynamic_covariate, id, landmark, n_pts)
      } else {
        NULL
      }
    },
    error = function(e) {
      warning(
        "Could not compute a model-based trajectory curve for plotting: ",
        conditionMessage(e)
      )
      NULL
    }
  )
}
