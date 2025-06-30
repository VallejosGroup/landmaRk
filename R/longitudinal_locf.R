#' Makes predictions from an LOCF process
#'
#' @param x A dataframe
#' @param newdata A data frame containing static covariates and individual
#'   IDs
#' @param subject Name of the column in newdata where individual IDs are stored.
#'
#' @returns A vector of predictions.
#'
#' @examples
predict_locf_ <- function(x, newdata, subject) {
  # pprob contains probability of observation belonging to a certain cluster
  # But it is possible that pprob does not contain predictions for some individuals
  # in newdata. That is because those individuals have not been use in training.
  # We will augment pprob using the sample average for individuals not used in
  # training first, find out the largest cluster.
  pprob <- x$pprob
  mode_cluster <- as.integer(names(sort(-table(pprob$class)))[1])
  # Allocation of clusters for prediction
  if (nrow(newdata) == nrow(pprob)) {
    cluster_allocation <- data.frame(
      id = pprob[, subject],
      cluster = pprob$class
    )
  } else {
    # Create a dataframe for individuals not used in training with the most common cluster
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
      nrow(predictions$pred),
      nrow(newdata)
    ))
  }
  # Choose correct cluster for prediction
  if (avg == FALSE) {
    #### predictions <- predictions$pred * model.matrix(~as.factor(cluster)-1,data=cluster_allocation)
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
  # Store predictions in Landmarking object
  names(predictions) <- newdata[, subject]
  predictions
}
