#' Concordance index for competing risks
#'
#' Assess discriminative performance of predictions obtained from a conventional
#' or competing risks time-to-event model using time-dependent concordance
#' index.
#'
#' Uses the proportion of correctly ordered risk pairs for the event
#' \eqn{k}, based on the predicted risk of the event up
#' to time \eqn{\tau}.
#'
#' \deqn{C_k(\tau) = \frac{\sum_{i=1}^N \sum_{j=1}^N (A_{ij} + B_{ij}) \cdot Q_{ij} \cdot N_i^k(\tau)}{\sum_{i=1}^N \sum_{j=1}^N (A_{ij} + B_{ij}) \cdot N_i^k(\tau)}}
#'
#' A == risk ordering of patients, small time means patient 'i' at higher risk
#' than patient 'j' experiencing event of interest \eqn{A[i,j] = 0} for tied
#' event times.
#'
#' B == risk ordering of patients, large time for patient 'i' means lower risk
#' than patient 'j' if not experienced the event of interest. Ties are included
#' in B
#'
#' Q == the risk ordering of the subjects, i.e., is subject i assigned a higher
#' risk by the model than the subject j, for event \eqn{E_k} until time \eqn{t}.
#' \eqn{Q[i,j] = 0} for tied predictions.
#'
#' N_t == number of subjects with survival time < time point and experience
#' event of interest. Tied event times are included
#'
#'
#' @param predictions Numeric vector of model predictions.
#' @param time Numeric vector describing the time to the event of interest or
#'   censoring.
#' @param cens.code Value used to denote censoring in \code{status}. Defaults to
#'   0.
#' @param status Vector of censoring status.
#' @param tau Time c-index is evaluated.
#' @param cause Event of interest.
#' @param method \code{'survival'} if the predictions are survival probabilities
#'   or \code{'cifs'} if they are cumulative incidence functions
#'
#' @references Ahuja K, Schaar M van der. Joint Concordance Index. Published
#'   online August 17, 2019. \doi{10.48550/arXiv.1810.11207}
#'
#' @return Concordance index value.

.CIndexCRisks <- function(
  predictions,
  time,
  cens.code = 0,
  status,
  cause,
  tau,
  method = c("survival", "cifs")
) {
  method <- match.arg(method)

  censor.status <- ifelse(status == cens.code, 0, 1)

  if (method == "survival") {
    predictions <- 1 - predictions
  }

  n <- length(predictions)
  A <- matrix(0, nrow = n, ncol = n)
  B <- matrix(0, nrow = n, ncol = n)
  Q <- matrix(0, nrow = n, ncol = n)
  N_t <- matrix(0, nrow = n, ncol = n)
  Num_mat <- matrix(0, nrow = n, ncol = n)
  Den_mat <- matrix(0, nrow = n, ncol = n)

  Num <- 0
  Den <- 0
  for (i in 1:n) {
    A[i, which(time[i] < time)] <- 1
    B[
      i,
      intersect(
        intersect(
          which((time[i] >= time)),
          which(status != cause)
        ),
        which(censor.status == 1)
      )
    ] <- 1
    Q[i, which(predictions[i] > predictions)] <- 1
  }

  for (i in 1:n) {
    if (time[i] <= tau && status[i] == cause && censor.status[i] == 1) {
      N_t[i, ] <- 1
    }
  }

  Num_mat <- (A + B) * Q * N_t
  Den_mat <- (A + B) * N_t

  Num <- sum(Num_mat)
  Den <- sum(Den_mat)

  Num / Den
}

#' Binary Brier Score
#'
#' Computes the Binary Brier Score (BBS) for binary outcomes at a specified
#' time point, \eqn{\tau}.
#'
#' The BBS is a measure of the accuracy of probabilistic predictions for binary
#' outcomes, where the predictions are the predicted probabilities of the event
#' of interest occurring by time \eqn{\tau}. The BBS is defined as the mean
#' squared difference between the predicted probabilities and the true outcome.
#' @inheritParams .CIndexCRisks
#' @param tau Time Brier score is evaluated.

.BinaryBrierScore <- function(predictions, time, status, tau, cause) {
  y_true <- ((time <= tau) * (status == cause))
  BS <- mean((predictions - y_true)^2)
  BS
}

# Validate inputs for performance_metrics function
.validate_performance_inputs <- function(
  x,
  landmarks,
  horizons,
  c_index,
  brier,
  auc_t,
  train,
  h_times
) {
  error_str <- NULL

  if (!inherits(x, "LandmarkAnalysis")) {
    error_str <- c(
      error_str,
      "@x must be an object of class LandmarkAnalysis"
    )
  }
  if (!is.numeric(landmarks)) {
    error_str <- c(error_str, "@landmarks must be a vector of numeric values")
  }
  if (!is.numeric(horizons)) {
    error_str <- c(error_str, "@horizons must be a vector of numeric values")
  }
  if (!(is.logical(c_index) && length(c_index) == 1L && !is.na(c_index))) {
    error_str <- c(error_str, "@c_index must be a length-1, non-NA logical")
  }
  if (!(is.logical(brier) && length(brier) == 1L && !is.na(brier))) {
    error_str <- c(error_str, "@brier must be a length-1, non-NA logical")
  }
  if (!(is.logical(auc_t) && length(auc_t) == 1L && !is.na(auc_t))) {
    error_str <- c(error_str, "@auc_t must be a length-1, non-NA logical")
  }
  if (!(is.logical(train) && length(train) == 1L && !is.na(train))) {
    error_str <- c(error_str, "@train must be a length-1, non-NA logical")
  }
  if (length(landmarks) != length(horizons)) {
    error_str <- c(
      error_str,
      "@landmarks and @horizons must be of the same length"
    )
  }

  if (length(error_str) > 0) {
    stop(paste(error_str, collapse = ". "))
  }
}
