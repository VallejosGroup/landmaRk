# Design-matrix row(s) for `newdata`. Used to build contrast vectors for
# hazard ratios: for two rows differing only in one covariate, their
# difference x1 - x0 gives log(HR) = (x1-x0) %*% coef(fit) and
# Var(log(HR)) = (x1-x0) %*% vcov(fit) %*% (x1-x0), for any formula
# (including interactions/splines), without hard-coding its functional form.
.coxph_design_matrix <- function(fit, newdata) {
  # terms in the formula expect response
  tt <- stats::delete.response(stats::terms(fit))
  for (special in c("strata", "cluster")) {
    untangled <- survival::untangle.specials(tt, special)
    if (length(untangled$terms) > 0) {
      tt <- stats::drop.terms(tt, untangled$terms, keep.response = FALSE)
    }
  }
  X <- suppressWarnings(
    stats::model.matrix(tt, data = newdata, xlev = fit$xlevels)
  )
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  X[, names(stats::coef(fit)), drop = FALSE]
}

# log(HR) and its SE for a `contrast`-unit increase in `covariate`, with
# every other covariate in `reference` fixed. `reference` must supply a value
# for every variable `fit`'s formula needs other than `covariate` itself.
.hazard_ratio_contrast <- function(fit, covariate, contrast, reference) {
  row0 <- row1 <- reference
  row0[[covariate]] <- 0
  row1[[covariate]] <- contrast
  dx <- .coxph_design_matrix(fit, rbind(row0, row1))
  dx <- dx[2L, ] - dx[1L, ]
  beta <- stats::coef(fit)
  # A term collinear with strata(LM) (e.g. a pure f(LM) main effect)
  nonzero <- dx != 0
  list(
    log_hr = as.numeric(dx[nonzero] %*% beta[nonzero]),
    se = as.numeric(sqrt(dx %*% stats::vcov(fit) %*% dx))
  )
}

# Default reference values (mean if numeric, first level if a factor) for
# every variable `fit`'s formula needs other than `covariate`/`exclude`
# (`exclude` is typically "LM", which callers supply explicitly per
# landmark instead of defaulting). Values are taken from `data` (normally
# `data_static`); any variable not found there (e.g. a dynamic covariate,
# which is never a column of `data_static`) is instead looked up in
# `model_data`, the data frame the model was actually fitted on (e.g.
# `x@survival_datasets[[key]]` or `x@superlandmark_datasets[[model_name]]`).
.default_reference_row <- function(
  fit,
  covariate,
  exclude,
  data,
  model_data = NULL
) {
  vars <- setdiff(
    all.vars(stats::delete.response(stats::terms(fit))),
    c(covariate, exclude)
  )
  row <- lapply(vars, function(v) {
    col <- if (v %in% colnames(data)) {
      data[[v]]
    } else if (!is.null(model_data) && v %in% colnames(model_data)) {
      model_data[[v]]
    } else {
      return(NULL)
    }
    if (is.numeric(col)) {
      mean(col, na.rm = TRUE)
    } else {
      factor(levels(as.factor(col))[1], levels = levels(as.factor(col)))
    }
  })
  names(row) <- vars
  row <- row[!vapply(row, is.null, logical(1))]
  if (length(row) == 0) {
    return(data.frame(row.names = 1L))
  }
  as.data.frame(row)
}

# .cutLM() is copied from cutLM() in the archived 'dynpred' package
# (https://github.com/cran/dynpred/blob/master/R/cutLM.R), by Hein Putter
# licensed GPL (>= 2).

#' Create landmark data set
#'
#' Create landmark data set from original data, which can be either in wide
#' or long format, see details.
#'
#' For a given landmark time point \code{LM}, patients who have reached the
#' event of interest (outcome) or are censored before or at \code{LM} are
#' removed. Administrative censoring is applied at the time horizon.
#' Time-varying covariates are evaluated at the landmark time point
#' \code{LM}.
#'
#' @param data Data frame from which to construct landmark dataset.
#' @param outcome List with items \code{time} and \code{status}, containing
#'   character strings identifying the names of time and status variables,
#'   respectively, of the survival outcome.
#' @param LM Scalar, the value of the landmark time point.
#' @param horizon Scalar, the value of the horizon. Administrative censoring
#'   is applied at \code{horizon}.
#' @param covs List with items \code{fixed} and \code{varying}, containing
#'   character strings specifying column names in the data containing
#'   time-fixed and time-varying covariates, respectively.
#' @param format Character string specifying whether the original data are
#'   in wide (default) or in long format.
#' @param id Character string specifying the column name in \code{data}
#'   containing the subject id; only needed if \code{format = "long"}.
#' @param rtime Character string specifying the column name in \code{data}
#'   containing the (running) time variable associated with the
#'   time-varying covariate(s); only needed if \code{format = "long"}.
#' @param right Boolean (default \code{TRUE}), indicating if the intervals
#'   for the time-varying covariates are closed on the right (and open on
#'   the left) or vice versa, see \code{\link{cut}}.
#'
#' @returns A landmark data set, containing the outcome and the values of
#'   time-fixed and time-varying covariates taken at the landmark time
#'   point. The value of the landmark time point is stored in column
#'   \code{LM}.
#'
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @references van Houwelingen HC, Putter H (2012). Dynamic Prediction in
#'   Clinical Survival Analysis. Chapman & Hall.
#' @noRd
.cutLM <- function(
  data,
  outcome,
  LM,
  horizon,
  covs,
  format = c("wide", "long"),
  id,
  rtime,
  right = TRUE
) {
  format <- match.arg(format)
  if (format == "wide") {
    LMdata <- data
    if (!is.null(covs$varying)) {
      LMdata[[covs$varying]] <- 1 - as.numeric(LMdata[[covs$varying]] > LM)
    }
  } else {
    if (missing(id)) {
      stop("argument 'id' should be specified for long format data")
    }
    if (missing(rtime)) {
      stop("argument 'rtime' should be specified for long format data")
    }
    ord <- order(data[[id]], data[[rtime]])
    data <- data[ord, ]
    ids <- unique(data[[id]])
    n <- length(ids)
    # initialize LMdata; copy first row of each subject
    LMdata <- data[which(!duplicated(data[[id]])), ]
    for (i in 1:n) {
      wh <- which(data[[id]] == ids[i])
      di <- data[wh, ]
      idx <- cut(LM, c(data[[rtime]][wh], Inf), right = right, labels = FALSE)
      if (!is.na(idx)) {
        LMdata[i, ] <- di[idx, ]
      } else {
        LMdata[i, ] <- di[1, ]
        LMdata[[covs$varying]][i] <- NA
        LMdata[[rtime]][i] <- NA
      }
    }
  }
  LMdata <- LMdata[LMdata[[outcome$time]] > LM, ]
  if (format == "long") {
    LMdata <- LMdata[!is.na(LMdata[[id]]), ]
  }
  # apply administrative censoring at horizon
  LMdata[outcome$status] <- LMdata[[outcome$status]] *
    as.numeric(LMdata[[outcome$time]] <= horizon)
  LMdata[outcome$time] <- pmin(as.vector(LMdata[[outcome$time]]), horizon)
  LMdata$LM <- LM
  if (format == "long") {
    cols <- match(
      c(
        id,
        outcome$time,
        outcome$status,
        covs$fixed,
        covs$varying,
        rtime,
        "LM"
      ),
      names(LMdata)
    )
  } else {
    cols <- match(
      c(outcome$time, outcome$status, covs$fixed, covs$varying, "LM"),
      names(LMdata)
    )
  }
  LMdata[, cols]
}
