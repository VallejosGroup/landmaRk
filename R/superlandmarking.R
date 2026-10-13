# Labels superlandmarking model according to landmark times and horizon
.superlandmark_model_name <- function(landmarks, w) {
  paste0(paste(sort(landmarks), collapse = ","), "-", w)
}

# Cluster() term to be added to formula in coxph() call
# to perform superlandmarking analysis
.add_cluster_term <- function(formula, ids) {
  rhs <- as.character(formula)[3]
  if (!grepl("cluster\\(", rhs)) {
    formula <- stats::as.formula(
      paste(
        as.character(formula)[2],
        as.character(formula)[1],
        rhs,
        "+ cluster(",
        ids,
        ")"
      )
    )
  }
  formula
}

# creates superlandmark dataframe
.create_superlandmark_dataframe <- function(
  x,
  landmarks,
  w,
  dynamic_covariates,
  include_clusters,
  validation_fold,
  train
) {
  # Iterate over landmark times to call .create_survival_dataframe()
  landmark_dfs <- lapply(landmarks, function(landmark) {
    survival_df <- .create_survival_dataframe(
      x,
      landmark,
      landmark + w,
      dynamic_covariates,
      include_clusters,
      censor_at_horizon = FALSE,
      validation_fold,
      train = train
    )
    fixed_covs <- setdiff(colnames(survival_df), c("event_time", "event_status"))
    landmark_df <- .cutLM(
      data = survival_df,
      outcome = list(time = "event_time", status = "event_status"),
      LM = 0,
      horizon = w,
      covs = list(fixed = fixed_covs, varying = NULL),
      format = "wide"
    )
    landmark_df$LM <- landmark
    landmark_df
  })

  # Stack LM-specific dataframes into a single dataframe
  dplyr::bind_rows(landmark_dfs)
}

#' Fits a survival model pooled across several landmark times (superlandmarking)
#'
#' @details
#'  This function implements "superlandmarking": a single Cox model is fitted to a
#'  "superdataset" subjects are present as many times as
#'  \code{landmarks} they are at risk. The landmark time is available to \code{formula} as
#'  the \code{LM} column, so that covariate effects and/or the baseline
#'  hazard can be allowed to vary with landmark time (via an interaction term with
#'  covariates f(LM), or \code{strata(LM)}).
#'
#'  Because each subject contributes one row per landmark at which they
#'  are still at risk, a \code{cluster(@ids)} term is added to
#'  \code{formula} automatically (unless already present).
#'
#' @param x An object of class \code{\link{LandmarkAnalysis}}.
#' @param formula A formula to be used in survival submodel fitting. Must
#'   reference the \code{LM} column directly to model landmark-varying
#'   covariate effects (e.g. via an interaction with a function of
#'   \code{LM}) or a landmark-varying baseline hazard (via
#'   \code{strata(LM)}).
#' @param landmarks Numeric vector of at least two landmark times, for which
#'   risk sets must already have been computed via
#'   \code{\link{compute_risk_sets}}.
#' @param w Scalar prediction window. Each landmark time \eqn{s} in
#'   \code{landmarks} uses the (fixed) horizon \eqn{s+w}.
#' @param dynamic_covariates Vector of time-varying covariates to be used in
#'   the survival model.
#' @param include_clusters Boolean indicating whether to propagate cluster
#'   membership to survival analysis.
#' @param validation_fold If positive, cross-validation fold where model is
#'   fitted. If 0 (default), model fitting is performed on the complete
#'   dataset.
#'
#' @returns An object of class \code{\link{LandmarkAnalysis}}.
#' @export
#'
#' @examples
#' data(epileptic)
#' epileptic_dfs <- split_wide_df(
#'   epileptic,
#'   ids = "id", times = "time",
#'   static = c("with.time", "with.status", "treat", "age", "gender", "learn.dis"),
#'   dynamic = c("dose"),
#'   measurement_name = "value"
#' )
#' x <- LandmarkAnalysis(
#'   data_static = epileptic_dfs$df_static,
#'   data_dynamic = epileptic_dfs$df_dynamic,
#'   event_indicator = "with.status",
#'   ids = "id", event_time = "with.time",
#'   times = "time", measurements = "value"
#' ) |>
#'   compute_risk_sets(c(90, 180, 365.25)) |>
#'   fit_superlandmark(
#'     formula = survival::Surv(event_time, event_status) ~
#'       treat * splines::ns(LM, 2) + age + strata(LM),
#'     landmarks = c(90, 180, 365.25),
#'     w = 2 * 365.25
#'   )
setGeneric(
  "fit_superlandmark",
  function(
    x,
    formula,
    landmarks,
    w,
    dynamic_covariates = c(),
    include_clusters = FALSE,
    validation_fold = 0
  ) {
    standardGeneric("fit_superlandmark")
  }
)

#' Fits a survival model pooled across several landmark times (superlandmarking)
#'
#' @inheritParams fit_superlandmark
#'
#' @returns An object of class \code{\link{LandmarkAnalysis}}.
#' @export
#'
#' @examples
setMethod(
  "fit_superlandmark",
  "LandmarkAnalysis",
  function(
    x,
    formula,
    landmarks,
    w,
    dynamic_covariates = c(),
    include_clusters = FALSE,
    validation_fold = 0
  ) {
    if (!is.numeric(landmarks) || length(landmarks) < 2) {
      stop("@landmarks must be a numeric vector of at least two landmark times.")
    }
    if (!is.numeric(w) || length(w) != 1 || w <= 0) {
      stop("@w must be a single, positive numeric prediction window.")
    }
    for (landmark in landmarks) {
      .check_riskset_survival(x, landmark)
    }

    model_name <- .superlandmark_model_name(landmarks, w)

    superdataset <- .create_superlandmark_dataframe(
      x,
      landmarks,
      w,
      dynamic_covariates,
      include_clusters,
      validation_fold,
      train = TRUE
    )

    formula <- .add_cluster_term(formula, x@ids)

    x@superlandmark_datasets[[model_name]] <- superdataset
    x@superlandmark_fits[[model_name]] <- .fit_coxph_survival(
      formula,
      data = superdataset
    )

    x
  }
)

#' Make predictions from a superlandmark survival model
#'
#' @inheritParams fit_superlandmark
#'
#' @returns An object of class \code{\link{LandmarkAnalysis}}.
#' @export
#'
#' @examples
setGeneric(
  "predict_superlandmark",
  function(
    x,
    landmarks,
    w,
    dynamic_covariates = c(),
    include_clusters = FALSE,
    validation_fold = 0
  ) {
    standardGeneric("predict_superlandmark")
  }
)

#' Make predictions from a superlandmark survival model
#'
#' @inheritParams predict_superlandmark
#'
#' @returns An object of class \code{\link{LandmarkAnalysis}}.
#' @export
#'
#' @examples
setMethod(
  "predict_superlandmark",
  "LandmarkAnalysis",
  function(
    x,
    landmarks,
    w,
    dynamic_covariates = c(),
    include_clusters = FALSE,
    validation_fold = 0
  ) {
    model_name <- .superlandmark_model_name(landmarks, w)
    if (!(model_name %in% names(x@superlandmark_fits))) {
      stop(
        "No superlandmark model has been fitted for landmarks ",
        paste(landmarks, collapse = ", "),
        " and window ",
        w,
        ". Call fit_superlandmark() first."
      )
    }

    x@superlandmark_predictions[[model_name]] <- survival::survfit(
      x@superlandmark_fits[[model_name]],
      newdata = x@superlandmark_datasets[[model_name]]
    )

    # Out of sample predictios
    if (validation_fold > 0) {
      x@superlandmark_datasets_test[[model_name]] <-
        .create_superlandmark_dataframe(
          x,
          landmarks,
          w,
          dynamic_covariates,
          include_clusters,
          validation_fold,
          train = FALSE
        )
      x@superlandmark_predictions_test[[model_name]] <- survival::survfit(
        x@superlandmark_fits[[model_name]],
        newdata = x@superlandmark_datasets_test[[model_name]]
      )
    }

    x
  }
)

#' Plot the hazard ratio of a covariate as a function of landmark time
#'
#' @details
#'  For every landmark time, the hazard ratio associated with a
#'  \code{contrast}-unit increase in \code{covariate} is computed as
#'  \eqn{\log(HR) = (x_1 - x_0)^\top \hat\beta}.
#'
#'  Any other covariate present in the formula is held fixed at the mean
#'  (numeric) or first level (factor) of \code{x@@data_static}, unless
#'  overridden via \code{newdata}.
#'
#' @param x An object of class \code{\link{LandmarkAnalysis}}.
#' @param covariate Character name of the (numeric) covariate to compute the
#'   hazard ratio for.
#' @param landmarks Numeric vector of landmark times. For \code{type =
#'   "crude"}, a fit must already exist in \code{x@@survival_fits} for each
#'   landmark (via \code{\link{fit_survival}}) at window \code{w}; for
#'   \code{type = "supermodel"}, a superlandmark model must already exist
#'   for this \code{landmarks}/\code{w} pair (via
#'   \code{\link{fit_superlandmark}}).
#' @param w Scalar prediction window, as passed to \code{\link{fit_survival}}
#'   (horizons \code{landmarks + w}) and/or \code{\link{fit_superlandmark}}.
#' @param type Whether to plot the crude per-landmark estimates, the pooled
#'   supermodel curve, or both.
#' @param contrast Numeric increase in \code{covariate} to compute the
#'   hazard ratio for (defaults to 1, i.e. the usual per-unit hazard ratio).
#' @param newdata A single-row \code{data.frame} giving values for every
#'   other covariate in the model (excluding \code{covariate} and \code{LM}).
#'   Defaults to the mean/first level of each from \code{x@@data_static}.
#' @param level Confidence level for the pointwise intervals (default 0.95).
#'
#' @returns A \code{\link[ggplot2]{ggplot}} object.
#' @export
#'
#' @examples
setGeneric(
  "plot_hazard_ratio",
  function(
    x,
    covariate,
    landmarks,
    w,
    type = c("crude", "supermodel", "both"),
    contrast = 1,
    newdata = NULL,
    level = 0.95
  ) {
    standardGeneric("plot_hazard_ratio")
  }
)

#' Plot the hazard ratio of a covariate as a function of landmark time
#'
#' @inheritParams plot_hazard_ratio
#'
#' @returns A \code{\link[ggplot2]{ggplot}} object.
#' @export
#' @importFrom ggplot2 ggplot geom_step geom_line aes labs theme_minimal
#'
#' @examples
setMethod(
  "plot_hazard_ratio",
  "LandmarkAnalysis",
  function(
    x,
    covariate,
    landmarks,
    w,
    type = c("crude", "supermodel", "both"),
    contrast = 1,
    newdata = NULL,
    level = 0.95
  ) {
    type <- match.arg(type)
    z <- stats::qnorm(1 - (1 - level) / 2)
    landmarks <- sort(landmarks)

    crude <- NULL
    if (type %in% c("crude", "both")) {
      crude <- lapply(landmarks, function(s) {
        key <- paste0(s, "-", s + w)
        fit <- x@survival_fits[[key]]
        if (is.null(fit)) {
          stop(
            "No crude survival model has been fitted for landmark ", s,
            " and window ", w, ". Call fit_survival() first."
          )
        }
        reference <- if (is.null(newdata)) {
          .default_reference_row(fit, covariate, character(0), x@data_static)
        } else {
          newdata
        }
        hr <- .hazard_ratio_contrast(fit, covariate, contrast, reference)
        data.frame(s = s, coef = hr$log_hr, se = hr$se)
      }) |>
        dplyr::bind_rows()
    }

    supermodel <- NULL
    if (type %in% c("supermodel", "both")) {
      model_name <- .superlandmark_model_name(landmarks, w)
      fit <- x@superlandmark_fits[[model_name]]
      if (is.null(fit)) {
        stop(
          "No superlandmark model has been fitted for landmarks ",
          paste(landmarks, collapse = ", "), " and window ", w,
          ". Call fit_superlandmark() first."
        )
      }
      reference <- if (is.null(newdata)) {
        .default_reference_row(fit, covariate, "LM", x@data_static)
      } else {
        newdata
      }
      fine_grid <- seq(min(landmarks), max(landmarks), length.out = 200)
      supermodel <- lapply(fine_grid, function(s) {
        reference$LM <- s
        hr <- .hazard_ratio_contrast(fit, covariate, contrast, reference)
        data.frame(s = s, coef = hr$log_hr, se = hr$se)
      }) |>
        dplyr::bind_rows()
    }

    p <- ggplot2::ggplot()
    if (!is.null(crude)) {
      p <- p +
        ggplot2::geom_step(
          data = crude,
          ggplot2::aes(.data[["s"]], .data[["coef"]] + z * .data[["se"]]),
          linetype = "solid", linewidth = 0.3
        ) +
        ggplot2::geom_step(
          data = crude,
          ggplot2::aes(.data[["s"]], .data[["coef"]] - z * .data[["se"]]),
          linetype = "solid", linewidth = 0.3
        ) +
        ggplot2::geom_step(
          data = crude,
          ggplot2::aes(.data[["s"]], .data[["coef"]], linetype = "Crude"),
          linewidth = 0.8
        )
    }
    if (!is.null(supermodel)) {
      p <- p +
        ggplot2::geom_line(
          data = supermodel,
          ggplot2::aes(.data[["s"]], .data[["coef"]] + z * .data[["se"]]),
          linetype = "dashed", linewidth = 0.3
        ) +
        ggplot2::geom_line(
          data = supermodel,
          ggplot2::aes(.data[["s"]], .data[["coef"]] - z * .data[["se"]]),
          linetype = "dashed", linewidth = 0.3
        ) +
        ggplot2::geom_line(
          data = supermodel,
          ggplot2::aes(.data[["s"]], .data[["coef"]], linetype = "Supermodel"),
          linewidth = 0.8
        )
    }
    p +
      ggplot2::labs(
        x = "Landmark time",
        y = paste0("Log hazard ratio (", covariate, ")"),
        linetype = NULL
      ) +
      ggplot2::theme_minimal()
  }
)
