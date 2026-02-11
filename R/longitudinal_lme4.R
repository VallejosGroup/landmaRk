#' Makes predictions from an LME model fitted using \code{\link[lme4]{lmer}}
#'
#' @param x An object of class \code{\link[lme4]{merMod}}.
#' @param newdata A data frame containing static covariates and individual
#'   IDs
#' @param subject Name of the column in newdata where individual IDs are stored.
#' @param test Logical indicating whether to make predictions for the test set
#'   (make out of sample predictions). Defaults to FALSE
#' @param newdata_long A data frame containing longitudinal measurements for
#'   prediction.
#' @param ... Optional additional parameters
#'
#' @returns A vector of predictions
#'
#' @noRd
#'
#' @examples
.predict_lme4 <- function(
  x,
  newdata,
  subject,
  test = FALSE,
  newdata_long,
  ...
) {
  # In the test set, lme4_predictions does not include the random-effect
  # component, so we calculate it and add it to lme4_predictions.
  if (!test) {
    lme4_predictions <- predict(
      x,
      newdata = newdata,
      allow.new.levels = TRUE
    )
    return(lme4_predictions)
    # Outwith the training set, lme4 predictions are marginal predictions
    # We make BLUP predictions and add them to the marginal predictions
  } else {
    # First, calculate predictions using lme4::predict.merMod
    lme4_predictions <- predict(
      x,
      newdata = newdata,
      allow.new.levels = TRUE
    )
    sigma2 <- sigma(x)^2
    beta_hat <- lme4::fixef(x)
    X_newdata <- model.matrix(formula(x, fixed.only = TRUE), newdata_long)
    attributes(X_newdata)$assign <- NULL
    attributes(X_newdata)$contrasts <- NULL

    # Random-effect variance-covariance matrix
    Sigma_b <- as.matrix(unclass(lme4::VarCorr(x)[[1]]))
    attributes(Sigma_b)$correlation <- NULL
    attributes(Sigma_b)$stddev <- NULL
    Sigma_b <- kronecker(Matrix::Diagonal(nrow(newdata)), Sigma_b)

    # Random-effect design matrix
    re_terms <- lme4::findbars(formula(x))[[1]]
    Z_block <- model.matrix(
      reformulate(deparse(re_terms[[2]])),
      newdata_long
    )

    patient_index <- as.integer(factor(newdata_long[, subject]))
    blocks <- split(as.data.frame(Z_block), patient_index)
    block_matrices <- lapply(blocks, as.matrix)
    Z_newdata <- Matrix::bdiag(block_matrices)

    b_hat <- solve(crossprod(Z_newdata) / sigma2 + solve(Sigma_b)) %*%
      (Matrix::t(Z_newdata) / sigma2) %*%
      (newdata_long[, as.character(formula(x, fixed.only = TRUE)[[2]])] -
        as.vector(X_newdata %*% beta_hat))

    Z_block <- model.matrix(
      reformulate(deparse(re_terms[[2]])),
      newdata
    )
    b_hat <- matrix(as.vector(b_hat), ncol = ncol(Z_block), byrow = TRUE)
    return(lme4_predictions + rowSums(Z_block * b_hat))
  }
}
