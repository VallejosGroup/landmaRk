#' Makes predictions from an LME model fitted using \code{\link[lme4]{lmer}}
#'
#' @param x An object of class \code{\link[lme4]{merMod}}.
#' @param newdata A data frame containing static covariates and individual
#'   IDs
#' @param subject Name of the column in newdata where individual IDs are stored.
#' @param test Logical indicating whether to make predictions for the test set
#'   (make out of sample predictions). Defaults to FALSE
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
  ...
) {
  # First, calculate predictions using lme4::predict.merMod
  lme4_predictions <- predict(
    x,
    newdata = newdata,
    allow.new.levels = TRUE
  )
  # In the
  if (!test) {
    return(lme4_predictions)
    # Outwith the training set, lme4 predictions are marginal predictions
    # We make BLUP predictions and add them to the marginal predictions
  } else {
    sigma2 <- sigma(x)^2
    beta_hat <- lme4::fixef(x)
    X_newdata <- model.matrix(formula(x, fixed.only = TRUE), newdata)
    attributes(X_newdata)$assign <- NULL
    attributes(X_newdata)$contrasts <- NULL

    # Random-effect variance-covariance matrix
    Sigma_b <- as.matrix(unclass(lme4::VarCorr(x)[[1]]))
    attributes(Sigma_b)$correlation <- NULL
    attributes(Sigma_b)$stddev <- NULL

    # Random-effect design matrix
    re_terms <- lme4::findbars(formula(x))[[1]]
    Z_newdata <- model.matrix(
      reformulate(deparse(re_terms[[2]])),
      newdata
    )

    b_hat <- matrix(0, nrow = nrow(newdata), ncol = ncol(Z_newdata))
    for (i in 1:nrow(newdata)) {
      Z_i <- matrix(Z_newdata[i, ], nrow = 1)
      b_hat[i, ] <- solve(
        crossprod(Z_i) / sigma2 + solve(Sigma_b),
        (t(Z_i) / sigma2) %*%
          (newdata[i, as.character(formula(x, fixed.only = TRUE)[[2]])] -
            X_newdata[i, ] %*% beta_hat)
      )
    }

    return(lme4_predictions + rowSums(Z_newdata * b_hat))
  }
}
