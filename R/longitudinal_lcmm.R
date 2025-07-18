#' Fits an LCMM model
#'
#' @param formula Two-sided linear formula for the fixed effects in the LCMM.
#' @param data Data frame with data
#' @param mixture One-sided formula specifying the class-specific fixed effects.
#' @param random One-sided formula specifying the random effects.
#' @param subject Name of the column indicating individual ids in data
#' @param ng Number of clusters in the LCMM model
#' @param ... Additional arguments passed to the \code{\link[lcmm]{hlme}}
#'   function.
#' @seealso  [lcmm::hlme()]
#' @noRd
#' @returns An object of class hlme
#'
#' @examples
.fit_lcmm <- function(formula, data, mixture, random, subject, ng, ...) {
  model_init <- lcmm::hlme(
    formula,
    data = data,
    random = random,
    subject = subject,
    ng = 1,
    ...
  )
  model_fit <- lcmm::hlme(
    formula,
    data = data,
    mixture = mixture,
    random = random,
    subject = subject,
    ng = ng,
    B = model_init,
    ...
  )

  .check_lcmm_convergence(model_fit)

  model_fit$call$fixed <- formula
  model_fit$call$mixture <- mixture
  model_fit$call$random <- random

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
#' @param avg Boolean indicating whether to make predictions based on the
#'   most likely cluster (FALSE, default) or averaging over clusters (TRUE).
#' @param include_clusters Boolean indicating whether to include
#'   predicted class allocation in the predictions.
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
  include_clusters = FALSE
) {
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
  if (nrow(newdata) != nrow(pprob)) {
    warning(
      "Individuals ",
      paste(setdiff(newdata[, subject], pprob[, subject]), collapse = ", "),
      ", have not been used in LCMM model fitting. ",
      "Imputing values for those individuals"
    )
    # Assign individuals not included in model fitting to the biggest cluster
    pprob.extra <- data.frame(
      id = setdiff(newdata[, subject], pprob[, subject]),
      cluster = mode_cluster
    )

    # Compute the column means for the probability matrix (excluding id and class columns)
    prob_means <- colMeans(pprob[, -c(1, 2)])

    # Convert the column means into a dataframe and transpose it
    prob_means_df <- t(as.data.frame(prob_means))

    # Repeat the column means for each individual in pprob.extra
    repeated_means <- apply(prob_means_df, 2, rep, each = nrow(pprob.extra))

    # Combine the repeated means with pprob.extra
    pprob.extra <- cbind(pprob.extra, repeated_means)

    # Reset row names and column names to match the original pprob structure
    rownames(pprob.extra) <- NULL
    colnames(pprob.extra) <- colnames(pprob)
    pprob <- rbind(pprob, pprob.extra) |> arrange(get(subject))
  }

  predictions <- lcmm::predictY(
    x,
    newdata = newdata,
    var.time = var.time,
    marg = TRUE
  )
  if (nrow(predictions$pred) != nrow(newdata)) {
    stop(sprintf(
      paste(
        "lcmm::predictY produced %d predictions but expected %d predictions.\n",
        "Probable reason: static covariates contain missing data.\n"
      ),
      nrow(predictions$pred),
      nrow(newdata)
    ))
  }
  # Choose correct cluster for prediction
  if (avg == FALSE) {
    predictions <- rowSums(
      predictions$pred *
        model.matrix(
          ~ as.factor(pprob$class) - 1,
          data = as.data.frame(pprob$class)
        )
    )
  } else {
    predictions <- rowSums(predictions$pred * as.matrix(pprob[, -c(1, 2)]))
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
