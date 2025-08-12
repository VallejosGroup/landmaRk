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
    returndata = TRUE,
    ...
  )

  .check_lcmm_convergence(model_fit)

  model_fit$call$fixed <- formula
  model_fit$call$mixture <- mixture
  model_fit$call$random <- random
  model_fit$call$subject <- subject
  model_fit$call$ng <- ng
  model_fit$call$data <- data
  model_fit$call$B <- model_init

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
  include_clusters = FALSE,
  .k = 0,
  test = FALSE
) {
  hlme <- NULL
  # Step 1. we make predictions for individuals in the training set.

  # Step 1a. We estimate the random effects for individuals in the training set
  x$call[[1]] <- expr(hlme)
  predRE <- lcmm::predictRE(x, x$data, subject = subject, classpredRE = TRUE)

  # Step 1b. Find ids of individuals in the training set
  in_train_set <- unique(x$data[, subject])
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

  # Step 1c. Find class-specific predictions for individuals in the training set.
  if (!test) {
    predictions_step1 <- t(sapply(
      in_train_set,
      function(individual) {
        return(
          lcmm::predictY(
            x,
            newdata = newdata |> filter(get(subject) == individual),
            predRE = predRE |> filter(get(subject) == individual)
          )$pred
        )
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

  # Step 1d. Estimate most likely cluster, and cluster allocation probabilities,
  # for individuals in the training set
  predicted_class_step1 <- lcmm::predictClass(x, x$data)

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
    predictions <- rowSums(
      as.matrix(predictions[, -1]) *
        as.matrix(pprob[, -c(1, 2)])
    )
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
            ~ as.factor(pprob$class) - 1,
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
