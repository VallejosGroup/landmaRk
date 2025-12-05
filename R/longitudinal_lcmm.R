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
  mixture,
  random,
  subject,
  ng,
  rep = 1,
  classmb = ~1,
  maxiter = 500,
  ...
) {
  model_init <- lcmm::hlme(
    formula,
    data = data,
    random = random,
    subject = subject,
    ng = 1,
    maxiter = maxiter
  )

  hlme <- NULL
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
      returndata = TRUE,
      maxiter = maxiter,
      ...
    )
  } else {
    model_fit <- lcmm::gridsearch(
      m = hlme(
        formula,
        data = data,
        mixture = mixture,
        random = random,
        subject = subject,
        ng = ng,
        classmb = classmb,
        returndata = TRUE,
        maxiter = maxiter,
        ...
      ),
      rep = rep,
      maxiter = 1, # This argument is ignored by lcmm::gridsearch
      minit = model_init
    )
  }

  .check_lcmm_convergence(model_fit)

  model_fit$call$fixed <- formula
  model_fit$call$mixture <- mixture
  model_fit$call$random <- random
  model_fit$call$subject <- subject
  model_fit$call$ng <- ng
  model_fit$call$data <- data
  model_fit$call$B <- model_init
  model_fit$call$classmb <- classmb
  model_fit
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
#' @param newdata_long A data frame containing longitudinal measurements for prediction. Required when \code{test = TRUE} and either \code{avg = TRUE} or \code{include_clusters = TRUE}. Should include columns for subject IDs, time (\code{var.time}), and any time-varying covariates used in the model. Defaults to \code{NULL}.
#'
#' @returns If \code{include_clusters == FALSE}, a vector of predictions. If
#'   \code{include_clusters == TRUE}, a vector whose first column includes
#'   predictions and second column includes predicted class allocation
#'
#' @noRd
#'
#' @examples
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
  # Step 1. we make predictions for individuals in the training set.

  # Step 1a. We estimate the random effects for individuals in the training set
  x$call[[1]] <- expr(hlme)
  # Step 1b. Find ids of individuals in the training set
  in_train_set <- unique(x$data[, subject])
  if (!test) {
    predRE <- lcmm::predictRE(x, x$data, subject = subject, classpredRE = TRUE)

    if (length(unique(predRE[, subject])) != length(in_train_set)) {
      stop(sprintf(
        paste(
          "lcmm::predictRE produced %d predictions but expected %d predictions.\n",
          "Probable reason: static covariates contain missing data.\n"
        ),
        length(unique(predRE[, subject])),
        length(in_train_set)
      ))
    }
  }

  # Step 1c. Find class-specific predictions for individuals in the training set.
  if (!test) {
    predictions_step1 <- t(sapply(
      in_train_set,
      function(individual) {
        lcmm::predictY(
          x,
          newdata = newdata |> filter(get(subject) == individual),
          predRE = predRE |> filter(get(subject) == individual)
        )$pred
      }
    ))

    predictions_step1 <- as.data.frame(predictions_step1)
    predictions_step1[, subject] <- in_train_set
    predictions_step1 <- predictions_step1 |> relocate(subject)
    colnames(predictions_step1) <- c(
      subject,
      paste0("Ypred_class", 1:(ncol(predictions_step1) - 1))
    )
  }

  # Step 2. we make predictions for individuals outwith the training set.

  # Step 2a. Find ids of individuals outwith the training set
  not_in_train_set <- setdiff(
    unique(newdata[, subject]),
    in_train_set
  )

  # Step 2b. Find class-specific predictions for individuals outwith the training set.
  if (length(not_in_train_set) > 0) {
    predictions_step2 <- lcmm::predictY(
      x,
      newdata = newdata |>
        filter(!(get(subject) %in% in_train_set))
    )$pred
    predictions_step2 <- as.data.frame(predictions_step2)
    predictions_step2[, subject] <- not_in_train_set
    predictions_step2 <- predictions_step2 |> relocate(subject)
  }

  # pprob contains probabilities for subjects belonging to each certain cluster,
  # However posterior probabilities are unavailable for  individuals not
  # included in the model fitting.
  # We augment pprob using the sample average for individuals not used in
  # model fitting.
  pprob <- x$pprob
  # Find the largest cluster
  mode_cluster <- as.integer(names(sort(-table(pprob$class)))[1])
  # If there are individuals in newdata that had not been used in model fitting,
  # we augment pprob imputing the sample average in those individuals
  if (!test && (nrow(newdata) != nrow(pprob))) {
    warning(
      "Individuals ",
      paste(setdiff(newdata[, subject], pprob[, subject]), collapse = ", "),
      ", have not been used in LCMM model fitting. ",
      "Imputing values for those individuals"
    )
    # Assign individuals not included in model fitting to the biggest cluster
    pprob_extra <- data.frame(
      id = setdiff(newdata[, subject], pprob[, subject]),
      cluster = mode_cluster
    )

    # Compute the column means for the probability matrix
    # (excluding id and class columns)
    prob_means <- colMeans(pprob[, -c(1, 2)])

    # Convert the column means into a dataframe and transpose it
    prob_means_df <- t(as.data.frame(prob_means))

    # Repeat the column means for each individual in pprob_extra
    repeated_means <- apply(prob_means_df, 2, rep, each = nrow(pprob_extra))

    # Combine the repeated means with pprob_extra
    pprob_extra <- cbind(pprob_extra, repeated_means)

    # Reset row names and column names to match the original pprob structure
    rownames(pprob_extra) <- NULL
    colnames(pprob_extra) <- colnames(pprob)

    pprob <- rbind(pprob, pprob_extra) |> arrange(get(subject))
  } else if (test && include_clusters) {
    # In the test set, use lcmm::predictClass to estimate cluster allocation
    pprob <- lcmm::predictClass(x, newdata = newdata_long)
  }

  if (length(not_in_train_set) > 0) {
    if (test) {
      predictions <- predictions_step2 |>
        arrange(get(subject))
    } else {
      predictions <- rbind(
        predictions_step1,
        predictions_step2
      ) |>
        arrange(get(subject))
    }
  } else {
    predictions <- predictions_step1
  }
  # If avg == TRUE, we return an average weighted according to cluster
  # probabilities. If avg == FALSE, we return the prediction according to the
  # most likely cluster
  if (avg) {
    if (test) {
      class_predictions <- lcmm::predictClass(
        x,
        newdata_long,
        subject = subject
      )
      predictions <- rowSums(class_predictions[, -c(1, 2)] * predictions[, -1])
      names(predictions) <- NULL
    } else {
      predictions <- rowSums(
        as.matrix(predictions[, -1]) *
          as.matrix(pprob[, -c(1, 2)])
      )
    }
  } else {
    if (test) {
      model_matrix_aux <- pprob |>
        inner_join(newdata, by = subject) |>
        select(starts_with("prob")) |>
        as.matrix()
      predictions <- rowSums(as.matrix(predictions[, -1]) * model_matrix_aux)
    } else {
      predictions <- rowSums(
        as.matrix(predictions[, -1]) *
          model.matrix(
            ~ factor(pprob$class, levels = as.character(1:x$ng)) - 1,
            data = as.data.frame(pprob$class)
          )
      )
    }
  }

  # Store predictions in LandmarkAnalysis object
  names(predictions) <- newdata[, subject]

  if (include_clusters) {
    # Append class labels
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
#' @export
#'
#' @examples
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
