#' Fits an LCMM model
#'
#' @param formula Two-sided linear formula for the fixed effects in the LCMM.
#' @param data Data frame with data
#' @param mixture One-sided formula specifying the class-specific fixed effects.
#' @param subject Name of the column indicating individual ids in data
#' @param ng Number of clusters in the LCMM model
#' @param ... Additional arguments passed to the \code{\link[lcmm]{hlme}}
#'   function.
#' @seealso  [lcmm::hlme()]
#'
#' @returns An object of class hlme
#'
#' @examples
fit_lcmm_ <- function(formula, data, mixture, subject, ng, ...) {
  model_init <- lcmm::hlme(formula, data = data, subject = subject, ng = 1)
  model_fit <- lcmm::hlme(formula,
    data = data,
    mixture = mixture,
    subject = subject,
    ng = ng,
    B = model_init,
    ...
  )

  model_fit$call$fixed <- formula
  model_fit$call$mixture <- mixture

  model_fit
}

#' Makes predictions from an LCMM model
#'
#' @param x An object of class \code{\link[lcmm]{hlme}}.
#' @param newdata A data frame containing static covariates and individual
#'   IDs
#' @param subject Name of the column in newdata where individual IDs are stored.
#' @param avg Boolean indicating whether to make predictions based on the
#'   most likely cluster (FALSE, default) or averaging over clusters (TRUE).
#'
#' @returns A vector of predictions.
#'
#' @examples
predict_lcmm_ <- function(x, newdata, subject, avg = FALSE) {
  # pprob contains probability of observation belonging to a certain cluster
  # But it is possible that pprob does not contain predictions for some individuals
  # in newdata. That is because those individuals have not been use in training.
  # We will augment pprob using the sample average for individuals not used in
  # training first, find out the largest cluster.
  pprob <- x$pprob
  mode_cluster <- as.integer(names(sort(-table(pprob$class)))[1])
  # Allocation of clusters for prediction
  if (nrow(newdata) == nrow(pprob)) {
    cluster_allocation <- data.frame(id = pprob[, subject], cluster = pprob$class)
  } else {
    # Create a dataframe for individuals not used in training with the most common cluster
    pprob.extra <- data.frame(id = setdiff(newdata[, subject], pprob[, subject]), cluster = mode_cluster)

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
    ### cluster_allocation <- rbind(
    ###   data.frame(id = pprob[,subject], cluster = pprob$class),
    ###   data.frame(id = setdiff(newdata$patient.id, pprob[,subject]), cluster = mode_cluster)
    ### ) |> arrange(id) |> select(-id)
  }
  ### rownames(cluster_allocation) <- newdata$patient.id
  # Make predictions with lcmm package
  predictions <- lcmm::predictY(x, newdata = newdata)
  if (nrow(predictions$pred) != nrow(newdata)) {
    stop(sprintf(
      "lcmm::predictY produced %d predictions but expected %d. Probable reason: static covariates contain missing data.",
      nrow(predictions$pred), nrow(newdata)
    ))
  }
  # Choose correct cluster for prediction
  if (avg == FALSE) {
    #### predictions <- predictions$pred * model.matrix(~as.factor(cluster)-1,data=cluster_allocation)
    predictions <- rowSums(predictions$pred * model.matrix(~ as.factor(pprob$class) - 1, data = as.data.frame(pprob$class)))
  } else {
    predictions <- rowSums(predictions$pred * as.matrix(pprob[, -c(1, 2)]))
  }
  # Store predictions in Landmarking object
  names(predictions) <- newdata[, subject]
  predictions
}
