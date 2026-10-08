#' Fits an LCMM model
#'
#' @param formula Two-sided linear formula for the fixed effects in the LCMM.
#' @param data Data frame with data
#' @param mixture One-sided formula specifying the class-specific fixed effects.
#' @param random One-sided formula specifying the random effects.
#' @param subject Name of the column indicating individual ids in data
#' @param ng Number of clusters in the LCMM model.
#' @param rep Number of times the model fitting algorithm is run using grid
#'   search.
#' @param maxiter Maximum number of iterations for the LCMM optimiser.
#' @param cl Only used when \code{rep} > 1. Either an integer giving the
#'   number of cores, or a cluster created by
#'   \code{\link[parallel]{makeCluster}}, to parallelise
#'   \code{\link[lcmm]{gridsearch}}'s \code{rep} random restarts. This is a
#'   different axis of parallelism from \code{\link{fit_longitudinal}}'s
#'   \code{cores} argument, which parallelises across landmark times instead
#'   of across restarts within a single grid search; combining both
#'   is not supported. Requires the lcmm package to be attached
#'   (\code{library(lcmm)}). Defaults to \code{NULL} (sequential).
#' @param idiag Logical. If TRUE, the random effects covariance matrix is
#'   constrained to be diagonal. Passed to the initialization model and all
#'   subsequent model fits. Defaults to FALSE.
#' @param ... Additional arguments passed to the \code{\link[lcmm]{hlme}}
#'   function.
#' @seealso  [lcmm::hlme()]
#' @noRd
#' @returns An object of class hlme
#'
#' @examples
.fit_lcmm <- function(
  formula,
  data,
  mixture = NULL,
  random,
  subject,
  ng,
  rep = 1,
  classmb = ~1,
  maxiter = 500,
  cl = NULL,
  idiag = FALSE,
  ...
) {
  model_init <- lcmm::hlme(
    formula,
    data = data,
    random = random,
    subject = subject,
    ng = 1,
    idiag = idiag,
    returndata = TRUE,
    maxiter = maxiter
  )

  hlme <- NULL
  if (ng == 1) {
    model_fit <- model_init
  } else {
    if (rep == 1) {
      model_fit <- lcmm::hlme(
        formula,
        data = data,
        mixture = mixture,
        random = random,
        subject = subject,
        ng = ng,
        B = model_init,
        classmb = classmb,
        idiag = idiag,
        returndata = TRUE,
        maxiter = maxiter,
        ...
      )
    } else {
      # lcmm::gridsearch() re-evaluates its `m` call (rather than reusing a
      # value) both for each grid-search restart and for the final model, so
      # `hlme` must resolve as an unqualified function name; unlike the
      # lcmm::hlme() calls elsewhere in this file, that requires the lcmm
      # package to be attached, not just imported.
      if (!("package:lcmm" %in% search())) {
        stop(
          "Fitting an LCMM model with `rep` > 1 (grid search) requires the ",
          "lcmm package to be attached: please call library(lcmm) before ",
          "fit_longitudinal().\n"
        )
      }
      if (is.null(cl)) {
        model_fit <- lcmm::gridsearch(
          m = hlme(
            formula,
            data = data,
            mixture = mixture,
            random = random,
            subject = subject,
            ng = ng,
            classmb = classmb,
            idiag = idiag,
            returndata = TRUE,
            maxiter = 24000,
            ...
          ),
          rep = rep,
          maxiter = maxiter,
          minit = model_init
        )
      } else {
        model_fit <- .gridsearch_lcmm(
          formula = formula,
          data = data,
          mixture = mixture,
          random = random,
          subject = subject,
          ng = ng,
          classmb = classmb,
          idiag = idiag,
          model_init = model_init,
          rep = rep,
          maxiter = maxiter,
          cl = cl,
          ...
        )
      }
    }
  }

  .check_lcmm_convergence(model_fit)

  model_fit$call$fixed <- formula
  model_fit$call$mixture <- mixture
  model_fit$call$subject <- subject
  model_fit$call$ng <- ng
  model_fit$call$data <- data
  model_fit$call$B <- model_init
  if (ng > 1) {
    model_fit$call$classmb <- classmb
    model_fit$call$random <- random
  }
  model_fit
}

# Runs `rep` random restarts of lcmm::hlme() in parallel across a cluster,
# then re-fits from the best restart's estimates, reimplementing
# lcmm::gridsearch()'s `cl` branch.
.gridsearch_lcmm <- function(
  formula,
  data,
  mixture,
  random,
  subject,
  ng,
  classmb,
  idiag,
  model_init,
  rep,
  maxiter,
  cl,
  ...
) {
  if (model_init$conv != 1) {
    stop("Initial (minit) model did not converge; cannot run grid search.")
  }

  extra_args <- list(...)
  hlme_call <- function(fit_maxiter, B) {
    as.call(c(
      quote(lcmm::hlme),
      list(formula, data = data),
      list(
        mixture = mixture,
        random = random,
        subject = subject,
        ng = ng,
        classmb = classmb,
        idiag = idiag,
        returndata = TRUE,
        maxiter = fit_maxiter,
        B = B
      ),
      extra_args
    ))
  }

  owns_cluster <- !inherits(cl, "cluster")
  if (owns_cluster) {
    if (!is.numeric(cl)) {
      stop(
        "argument cl should be either a cluster or a numeric value ",
        "indicating the number of cores"
      )
    }
    cl <- parallel::makeCluster(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  parallel::clusterSetRNGStream(cl)

  # `random(model_init)` is never actually called: hlme() detects this
  # literal, unevaluated call passed as `B` and generates randomised
  # initial values from it internally.
  random_B <- as.call(list(quote(random), model_init))

  models <- parallel::parLapply(
    cl,
    seq_len(rep),
    function(i, hlme_call, random_B, maxiter) {
      eval(hlme_call(maxiter, random_B))
    },
    hlme_call = hlme_call,
    random_B = random_B,
    maxiter = maxiter
  )

  best <- models[[
    which.max(vapply(models, function(x) x$loglik, numeric(1)))
  ]]
  eval(hlme_call(24000, best$best))
}

.check_lcmm_convergence <- function(model_fit) {
  switch(
    as.character(model_fit$conv),
    "1" = message("LCMM model converged successfully."),
    "2" = warning("Maximum number of iterations reached without convergence."),
    "3" = message(
      "Convergence criteria satisfied with a partial Hessian matrix."
    ),
    stop("Problem occurred during optimisation of the LCMM model.")
  )
}

# Helper function to make class-specific predictions for a set of individuals
.make_class_predictions <- function(
  x,
  individuals,
  newdata,
  predRE,
  subject
) {
  predictions <- t(sapply(
    individuals,
    function(individual) {
      lcmm::predictY(
        x,
        newdata = newdata |> filter(get(subject) == individual),
        predRE = predRE |> filter(get(subject) == individual)
      )$pred
    }
  ))

  if (x$call$ng == 1) {
    predictions <- as.data.frame(t(predictions))
  } else {
    predictions <- as.data.frame(predictions)
  }
  predictions[, subject] <- individuals
  predictions <- predictions |> relocate(all_of(subject))
  predictions
}

#' Makes predictions from an LCMM model
#'
#' @param x An object of class \code{\link[lcmm]{hlme}}.
#' @param newdata A data frame containing static covariates and individual
#'   IDs
#' @param subject Name of the column in newdata where individual IDs are stored.
#' @param var.time Name of the column in newdata where time is recorded.
#' @param avg Logical indicating whether to make predictions based on the
#'   most likely cluster (FALSE, default) or averaging over clusters (TRUE).
#' @param include_clusters Logical indicating whether to include
#'   predicted class allocation in the predictions.
#' @param validation_fold If positive, cross-validation fold where model is
#'   fitted. If 0 (default), model fitting is performed in the complete dataset.
#' @param test Logical indicating whether to make predictions for the test set
#'   (make out of sample predictions). Defaults to FALSE
#' @param newdata_long A data frame containing longitudinal measurements for
#'   prediction. Required when \code{test = TRUE} and either \code{avg = TRUE}
#'   or \code{include_clusters = TRUE}. Should include columns for subject IDs,
#'   time (\code{var.time}), and any time-varying covariates used in the model.
#'   Defaults to \code{NULL}.
#'
#' @returns If \code{include_clusters == FALSE}, a vector of predictions. If
#'   \code{include_clusters == TRUE}, a vector whose first column includes
#'   predictions and second column includes predicted class allocation
#'
#' @noRd
#'
#' @examples
# Errors if the number of individuals with random effect predictions does not
# match the number expected.
.check_predRE_count <- function(predRE, subject, expected) {
  n_pred <- length(unique(predRE[, subject]))
  if (n_pred != expected) {
    stop(sprintf(
      paste(
        "lcmm::predictRE produced %d predictions but expected",
        "%d predictions.\n",
        "Probable reason: static covariates contain missing data.\n"
      ),
      n_pred,
      expected
    ))
  }
}

# Random effects predictions for individuals in the training set
.predict_re_train <- function(x, in_train_set, subject) {
  predRE <- lcmm::predictRE(
    x,
    x$data |> filter(get(subject) %in% in_train_set),
    subject = subject,
    classpredRE = TRUE
  )
  .check_predRE_count(predRE, subject, length(in_train_set))
  predRE
}

# Random effects predictions for the test set. Individuals without
# observations get zero random effects in every class.
.predict_re_test <- function(x, newdata, newdata_long, subject) {
  predRE <- lcmm::predictRE(
    x,
    newdata_long,
    subject = subject,
    classpredRE = TRUE
  )
  subjects_no_obs <- setdiff(
    newdata[, subject],
    unique(newdata_long[, subject])
  )
  if (length(subjects_no_obs) > 0) {
    re_cols <- setdiff(colnames(predRE), c(subject, "class"))
    predRE_default <- expand.grid(
      id = subjects_no_obs,
      class = 1:x$ng
    )
    colnames(predRE_default)[1] <- subject
    predRE_default[, re_cols] <- 0
    predRE <- rbind(predRE, predRE_default)
  }
  .check_predRE_count(predRE, subject, nrow(newdata))
  list(predRE = predRE, subjects_no_obs = subjects_no_obs)
}

# Converts the matrix of predictions made by sapply into a data frame with a
# leading subject column
.predictions_to_df <- function(predictions, ng, ids, subject) {
  if (ng == 1) {
    predictions <- as.data.frame(t(predictions))
  } else {
    predictions <- as.data.frame(predictions)
  }
  predictions[, subject] <- ids
  predictions |> relocate(all_of(subject))
}

# Class-specific predictions for individuals not used in model fitting, made
# without random effects predictions
.make_untrained_predictions <- function(x, individuals, newdata, subject) {
  predictions <- t(sapply(
    individuals,
    function(individual) {
      lcmm::predictY(
        x,
        newdata = newdata |> filter(get(subject) == individual)
      )$pred
    }
  ))
  .predictions_to_df(predictions, x$call$ng, individuals, subject)
}

# Class-specific predictions for all individuals in newdata
.make_all_class_predictions <- function(
  x,
  newdata,
  predRE,
  subject,
  in_train_set,
  not_in_train_set,
  test
) {
  if (test) {
    return(
      .make_class_predictions(x, not_in_train_set, newdata, predRE, subject) |>
        arrange(get(subject))
    )
  }
  predictions <- .make_class_predictions(
    x,
    in_train_set,
    newdata,
    predRE,
    subject
  )
  colnames(predictions) <- c(
    subject,
    paste0("Ypred_class", 1:(ncol(predictions) - 1))
  )
  if (length(not_in_train_set) == 0) {
    return(predictions)
  }
  predictions_step2 <- .make_untrained_predictions(
    x,
    not_in_train_set,
    newdata,
    subject
  )
  colnames(predictions_step2) <- colnames(predictions)
  rbind(predictions, predictions_step2) |>
    arrange(get(subject))
}

# Appends rows with sample-average class probabilities for individuals without
# observations
.append_default_class_probs <- function(
  probs,
  subjects_no_obs,
  mode_cluster
) {
  if (length(subjects_no_obs) == 0) {
    return(probs)
  }
  prob_means <- colMeans(probs[, -c(1, 2)])
  probs_default <- data.frame(
    id = subjects_no_obs,
    class = mode_cluster,
    matrix(
      rep(prob_means, each = length(subjects_no_obs)),
      nrow = length(subjects_no_obs),
      dimnames = list(NULL, names(prob_means))
    )
  )
  colnames(probs_default) <- colnames(probs)
  rbind(probs, probs_default)
}

# Augments pprob using the sample average for individuals not used in model
# fitting (posterior probabilities are unavailable for them), assigning them
# to the largest cluster.
.impute_pprob <- function(pprob, newdata, subject, mode_cluster) {
  missing_ids <- setdiff(newdata[, subject], pprob[, subject])
  warning(
    "Individuals ",
    paste(missing_ids, collapse = ", "),
    ", have not been used in LCMM model fitting. ",
    "Imputing values for those individuals"
  )
  pprob_extra <- data.frame(id = missing_ids, cluster = mode_cluster)

  # Column means of the probability matrix (excluding id and class columns),
  # repeated for each individual in pprob_extra
  prob_means_df <- t(as.data.frame(colMeans(pprob[, -c(1, 2)])))
  repeated_means <- apply(prob_means_df, 2, rep, each = nrow(pprob_extra))
  pprob_extra <- cbind(pprob_extra, repeated_means)

  rownames(pprob_extra) <- NULL
  colnames(pprob_extra) <- colnames(pprob)

  rbind(pprob, pprob_extra) |> arrange(get(subject))
}

# Posterior class probabilities for individuals in newdata
.get_pprob <- function(
  x,
  newdata,
  subject,
  test,
  include_clusters,
  newdata_long,
  subjects_no_obs,
  mode_cluster,
  avg = FALSE
) {
  pprob <- x$pprob |>
    filter(
      get(subject) %in%
        intersect(unique(newdata[, subject]), unique(x$data[, subject]))
    )
  if (!test && nrow(newdata) != nrow(pprob)) {
    return(.impute_pprob(pprob, newdata, subject, mode_cluster))
  }
  # Test subjects are absent from x$data, so x$pprob is empty for them. Class
  # allocation is needed whenever clusters are returned or a multi-class model
  # selects class-specific predictions (avg = FALSE).
  if (test && (include_clusters || (x$ng > 1 && !avg))) {
    # In the test set, use lcmm::predictClass to estimate cluster allocation
    pprob <- lcmm::predictClass(x, newdata = newdata_long, subject = subject)
    if (length(subjects_no_obs) > 0) {
      pprob <- .append_default_class_probs(
        pprob,
        subjects_no_obs,
        mode_cluster
      )
    }
    pprob <- pprob[match(newdata[, subject], pprob[, subject]), ]
  }
  pprob
}

# Weighted average of class-specific predictions using test-set class
# probabilities
.average_test_predictions <- function(
  x,
  predictions,
  newdata,
  newdata_long,
  subject,
  subjects_no_obs,
  mode_cluster
) {
  class_predictions <- lcmm::predictClass(x, newdata_long, subject = subject)
  class_predictions <- .append_default_class_probs(
    class_predictions,
    subjects_no_obs,
    mode_cluster
  )
  class_predictions <- class_predictions[
    match(newdata[, subject], class_predictions[, subject]),
  ]
  result <- rowSums(class_predictions[, -c(1, 2)] * predictions[, -1])
  names(result) <- NULL
  result
}

# Reduces class-specific predictions to a single prediction per individual
.reduce_class_predictions <- function(
  x,
  predictions,
  pprob,
  newdata,
  newdata_long,
  subject,
  avg,
  test,
  subjects_no_obs,
  mode_cluster
) {
  if (avg && test) {
    return(.average_test_predictions(
      x,
      predictions,
      newdata,
      newdata_long,
      subject,
      subjects_no_obs,
      mode_cluster
    ))
  }
  if (avg) {
    return(rowSums(
      as.matrix(predictions[, -1]) * as.matrix(pprob[, -c(1, 2)])
    ))
  }
  if (x$call$ng == 1) {
    return(predictions[, -1])
  }
  rowSums(
    as.matrix(predictions[, -1]) *
      model.matrix(
        ~ factor(pprob$class, levels = as.character(1:x$ng)) - 1,
        data = as.data.frame(pprob$class)
      )
  )
}

.predict_lcmm <- function(
  x,
  newdata,
  subject,
  var.time,
  avg = FALSE,
  include_clusters = FALSE,
  validation_fold = 0,
  classmb = NULL,
  test = FALSE,
  newdata_long = NULL
) {
  hlme <- NULL
  x$call[[1]] <- expr(hlme)

  in_train_set <- intersect(
    unique(newdata[, subject]),
    unique(x$data[, subject])
  )
  not_in_train_set <- setdiff(unique(newdata[, subject]), in_train_set)

  # Random effects
  subjects_no_obs <- NULL
  if (test) {
    re <- .predict_re_test(x, newdata, newdata_long, subject)
    predRE <- re$predRE
    subjects_no_obs <- re$subjects_no_obs
  } else {
    predRE <- .predict_re_train(x, in_train_set, subject)
  }

  predictions <- .make_all_class_predictions(
    x,
    newdata,
    predRE,
    subject,
    in_train_set,
    not_in_train_set,
    test
  )

  # Largest cluster
  mode_cluster <- as.integer(names(sort(-table(x$pprob$class)))[1])
  if (is.na(mode_cluster)) {
    mode_cluster <- 1L
  }
  pprob <- .get_pprob(
    x,
    newdata,
    subject,
    test,
    include_clusters,
    newdata_long,
    subjects_no_obs,
    mode_cluster,
    avg
  )

  # arrange()/rbind() above reorders rows by ascending subject id, which does
  # not generally match newdata's row order. Re-align predictions to
  # newdata's row order before names(predictions) <- newdata[, subject] is
  # assigned further down.
  predictions <- predictions[
    match(newdata[, subject], predictions[, subject]),
  ]

  # If avg == TRUE, we return an average weighted according to cluster
  # probabilities. If avg == FALSE, we return the prediction according to the
  # most likely cluster
  predictions <- .reduce_class_predictions(
    x,
    predictions,
    pprob,
    newdata,
    newdata_long,
    subject,
    avg,
    test,
    subjects_no_obs,
    mode_cluster
  )
  names(predictions) <- newdata[, subject]

  if (include_clusters) {
    predictions <- cbind(predictions, cluster = pprob[, "class"])
    predictions <- as.data.frame(predictions)
    predictions$cluster <- as.factor(predictions$cluster)
  }

  predictions
}

#' Checks convergence of lcmm models
#'
#' @param x An object of class \code{\link{LandmarkAnalysis}}.
#'
#' @return No return value. Issues a message if all lcmm models converged,
#' or a warning for each model that did not converge.
#'
#' @export
check_lcmm_convergence <- function(x) {
  if (!is(x, "LandmarkAnalysis")) {
    stop("x must be an object of class LandmarkAnalysis")
  } else if (length(x@longitudinal_fits) == 0) {
    stop(
      "Longitudinal submodels must be fitted before calling ",
      "check_lcmm_convergence()"
    )
  }
  num_models_not_converged <- 0
  for (landmark in names(x@longitudinal_fits)) {
    for (dynamic_covariate in names(x@longitudinal_fits[[landmark]])) {
      if (!(is(x@longitudinal_fits[[landmark]][[dynamic_covariate]], "hlme"))) {
        warning(paste0(
          "Longitudinal model for dynamic covariate ",
          dynamic_covariate,
          " at landmark time ",
          landmark,
          "was not fitted using LCMM."
        ))
      }
      conv_status <- x@longitudinal_fits[[landmark]][[dynamic_covariate]]$conv
      if (!(conv_status %in% c(1, 3))) {
        num_models_not_converged <- num_models_not_converged + 1
        msg <- paste0(
          "Model for dynamic covariate ",
          dynamic_covariate,
          " at landmark time ",
          landmark,
          " did not converge. ",
          switch(
            as.character(
              x@longitudinal_fits[[landmark]][[dynamic_covariate]]$conv
            ),
            "2" = "Maximum number of iterations were reached.",
            "Problem occured during optimisation."
          )
        )
        warning(msg)
      }
    }
  }
  if (num_models_not_converged == 0) {
    message("All longitudinal models converged.")
  }
}
