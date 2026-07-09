#' Plot longitudinal observations and predicted survival curve for one individual
#'
#' Produces a single, self-explanatory panel with a common time axis. To the
#' left of the landmark dashed line, the individual's observed longitudinal
#' measurements are shown, together with model-based context depending on the
#' longitudinal sub-model used at that landmark:
#' \itemize{
#'   \item \strong{LOCF} (or another summary measure): the last observed value
#'     is carried forward to the landmark (dashed segment).
#'   \item \strong{lme4}: the population-average trajectory (fixed effects
#'     only) and the individual's predicted trajectory (including their
#'     predicted random effects) are both drawn.
#'   \item \strong{lcmm}: the average trajectory of \emph{every} latent
#'     cluster is drawn, together with the individual's own predicted
#'     trajectory (each cluster's fixed effects plus the individual's
#'     predicted random effects, averaged across clusters using the
#'     individual's posterior class-membership probabilities). The
#'     individual's most likely cluster and their posterior probability of
#'     belonging to each cluster are noted in the plot title/subtitle.
#' }
#' In every case, the value that is actually fed into the survival sub-model
#' is highlighted at the landmark. To the right of the landmark dashed line,
#' the individual's predicted survival curve is shown on a secondary axis.
#' A legend identifies every series.
#'
#' @param x An object of class \code{\link{LandmarkAnalysis}}.
#' @param id Identifier of the individual to plot. Must match a value in the
#'   column \code{x@@ids}.
#' @param landmark Numeric landmark time.
#' @param dynamic_covariate Character name of the dynamic covariate to display.
#' @param horizon Numeric horizon time. If \code{NULL} (default), uses the
#'   single available horizon for \code{landmark}; errors when multiple
#'   horizons are available.
#' @param train Logical. If \code{TRUE} (default), uses in-sample predictions.
#'   If \code{FALSE}, uses out-of-sample predictions (requires
#'   \code{validation_fold > 0} in \code{\link{predict_survival}}).
#' @param ... Additional arguments (not currently used).
#'
#' @returns A \code{\link[ggplot2]{ggplot}} object.
#' @export
#' @importFrom ggplot2 ggplot geom_line geom_point geom_vline aes scale_y_continuous sec_axis xlab ggtitle theme_bw geom_segment scale_colour_manual
#'
#' @examples
setMethod(
  "plot",
  "LandmarkAnalysis",
  function(
    x,
    id,
    landmark,
    dynamic_covariate,
    horizon = NULL,
    train = TRUE,
    ...
  ) {
    .env <- NULL
    # ---- Input validation ----
    if (missing(id) || is.null(id)) {
      stop("@id is required.")
    }
    if (missing(landmark) || !is.numeric(landmark) || length(landmark) != 1) {
      stop("@landmark must be a single numeric value.")
    }
    if (!(landmark %in% x@landmarks)) {
      stop(paste(
        "Risk set for landmark time",
        landmark,
        "has not been computed."
      ))
    }
    if (missing(dynamic_covariate) || is.null(dynamic_covariate)) {
      stop("@dynamic_covariate is required.")
    }
    predictions_slot <- if (train) {
      x@survival_predictions
    } else {
      x@survival_predictions_test
    }
    if (length(predictions_slot) == 0) {
      stop(
        if (train) {
          "No survival predictions found. Call predict_survival() first."
        } else {
          "No out-of-sample survival predictions found. Call predict_survival() with validation_fold > 0."
        }
      )
    }

    # ---- Resolve horizon ----
    available <- grep(
      paste0("^", landmark, "-"),
      names(predictions_slot),
      value = TRUE
    )
    if (length(available) == 0) {
      stop(paste(
        "No survival predictions available for landmark time",
        landmark
      ))
    }
    if (is.null(horizon)) {
      if (length(available) == 1) {
        model_name <- available
        horizon <- as.numeric(sub(paste0(landmark, "-"), "", model_name))
      } else {
        stop(paste0(
          "Multiple horizons available for landmark ",
          landmark,
          ": ",
          paste(available, collapse = ", "),
          ". Please specify @horizon."
        ))
      }
    } else {
      model_name <- paste0(landmark, "-", horizon)
      if (!(model_name %in% names(predictions_slot))) {
        stop(paste(
          "No survival predictions for landmark",
          landmark,
          "and horizon",
          horizon
        ))
      }
    }

    # ---- Longitudinal observations ----
    obs_data <- x@data_dynamic[[dynamic_covariate]] |>
      dplyr::filter(.data[[x@ids]] == .env$id, .data[[x@times]] <= landmark)

    if (nrow(obs_data) == 0) {
      stop(paste(
        "No observations found for individual",
        id,
        "at or before landmark time",
        landmark
      ))
    }

    # ---- Predicted value at landmark (optional) ----
    long_preds <- if (train) {
      x@longitudinal_predictions[[as.character(landmark)]][[dynamic_covariate]]
    } else {
      x@longitudinal_predictions_test[[as.character(landmark)]][[
        dynamic_covariate
      ]]
    }
    pred_at_landmark <- if (is.null(long_preds)) {
      NULL
    } else if (is.data.frame(long_preds)) {
      as.numeric(long_preds[as.character(id), 1L])
    } else {
      as.numeric(long_preds[as.character(id)])
    }

    # ---- Model-based longitudinal trajectory (optional) ----
    # LOCF (and other summary measures) do not fit a model, so there is no
    # trajectory to compute; longitudinal_fits is left unpopulated in that
    # case (see fit_longitudinal()/predict_longitudinal()).
    long_fit <- x@longitudinal_fits[[as.character(landmark)]][[dynamic_covariate]]
    is_locf <- is.null(long_fit)
    trajectory <- if (is_locf) {
      NULL
    } else {
      .model_trajectory(long_fit, x, dynamic_covariate, id, landmark)
    }
    cluster <- trajectory$cluster
    probs <- trajectory$probs

    # For lcmm, every latent class gets its own average-trajectory label and
    # colour; for lme4 there is a single "population average" curve
    if (!is.null(trajectory)) {
      n_clusters <- nlevels(trajectory$avg_curves$cluster)
      cluster_labels <- if (!is.null(cluster)) {
        paste0("Cluster ", seq_len(n_clusters), " average trajectory")
      } else {
        "Population average trajectory"
      }
      trajectory$avg_curves$label <- cluster_labels[
        as.integer(trajectory$avg_curves$cluster)
      ]
      # Reuses hues from the validated categorical palette that are never
      # shown alongside a multi-cluster lcmm plot (LOCF's orange only
      # appears when there is no model fit at all)
      cluster_colours <- rep_len(
        c("#4a3aa7", "#eb6834", "#e87ba4", "#eda100"),
        n_clusters
      )
      names(cluster_colours) <- cluster_labels
    }

    # ---- Individual-specific survival curve ----
    dataset <- if (train) {
      x@survival_datasets[[model_name]]
    } else {
      x@survival_datasets_test[[model_name]]
    }
    sf <- predictions_slot[[model_name]]
    id_idx <- which(dataset[[x@ids]] == id)
    if (length(id_idx) == 0L) {
      stop(paste(
        "Individual",
        id,
        "not found in survival dataset for",
        model_name
      ))
    }

    n_pts <- 200L
    surv_times_rel <- seq(0, horizon - landmark, length.out = n_pts)
    sfit <- summary(sf, times = surv_times_rel)
    surv_mat <- matrix(sfit$surv, nrow = n_pts, ncol = nrow(dataset))
    surv_df <- data.frame(
      time = surv_times_rel + landmark,
      surv = surv_mat[, id_idx]
    )

    # ---- Y-axis range (longitudinal scale with a small buffer) ----
    y_vals <- c(
      obs_data[[x@measurements]],
      if (!is.null(pred_at_landmark)) pred_at_landmark,
      if (!is.null(trajectory)) {
        c(trajectory$avg_curves$average, trajectory$ind_curve$individual)
      }
    )
    y_raw <- range(y_vals, na.rm = TRUE)
    y_buf <- max(0.05 * diff(y_raw), 0.5)
    y_min <- y_raw[1L] - y_buf
    y_max <- y_raw[2L] + y_buf

    # Scale survival [0, 1] -> [y_min, y_max] for plotting on primary axis
    surv_df$surv_scaled <- surv_df$surv * (y_max - y_min) + y_min

    # ---- Colour legend, built up to match whichever series are shown ----
    palette <- c("Observed measurements" = "#2a78d6")
    if (!is.null(trajectory)) {
      palette <- c(
        palette,
        cluster_colours,
        "Individual predicted trajectory" = "#1baf7a"
      )
    }
    if (is_locf && !is.null(pred_at_landmark)) {
      palette <- c(palette, "LOCF carried forward" = "#eb6834")
    }
    if (!is.null(pred_at_landmark)) {
      palette <- c(palette, "Predicted value (landmark)" = "#e34948")
    }
    palette <- c(palette, "Survival probability" = "#008300")

    # ---- Assemble plot ----
    myplot <- ggplot2::ggplot() +
      # Survival curve
      ggplot2::geom_line(
        data = surv_df,
        ggplot2::aes(
          x = .data[["time"]],
          y = .data[["surv_scaled"]],
          colour = "Survival probability"
        ),
        linewidth = 0.8
      ) +
      # Observed longitudinal measurements
      ggplot2::geom_point(
        data = obs_data,
        ggplot2::aes(
          x = .data[[x@times]],
          y = .data[[x@measurements]],
          colour = "Observed measurements"
        )
      ) +
      # Vertical dashed line separating longitudinal from survival
      ggplot2::geom_vline(
        xintercept = landmark,
        linetype = "dashed",
        colour = "grey40"
      ) +
      # Model-based trajectory curves (LME population-average, or every
      # LCMM cluster-average, plus the individual-specific curve), when
      # available
      (if (!is.null(trajectory)) {
        list(
          ggplot2::geom_line(
            data = trajectory$avg_curves,
            ggplot2::aes(
              x = .data[["time"]],
              y = .data[["average"]],
              colour = .data[["label"]],
              group = .data[["cluster"]]
            ),
            linetype = "dashed",
            linewidth = 0.8
          ),
          ggplot2::geom_line(
            data = trajectory$ind_curve,
            ggplot2::aes(
              x = .data[["time"]],
              y = .data[["individual"]],
              colour = "Individual predicted trajectory"
            ),
            linewidth = 0.8
          )
        )
      }) +
      # Dashed segment from last observation to landmark summary value
      # (LOCF/summary measures only)
      (if (is_locf && !is.null(pred_at_landmark)) {
        last_obs_time <- obs_data |>
          dplyr::arrange(.data[[x@times]]) |>
          dplyr::slice_tail(n = 1L) |>
          dplyr::pull(x@times)
        ggplot2::geom_segment(
          data = data.frame(
            x = last_obs_time,
            xend = landmark,
            y = pred_at_landmark,
            yend = pred_at_landmark
          ),
          ggplot2::aes(
            x = .data[["x"]],
            xend = .data[["xend"]],
            y = .data[["y"]],
            yend = .data[["yend"]],
            colour = "LOCF carried forward"
          ),
          linetype = "dashed",
          linewidth = 0.8
        )
      }) +
      # Predicted value at landmark (diamond, highlighted; only if available)
      (if (!is.null(pred_at_landmark)) {
        ggplot2::geom_point(
          data = data.frame(t = landmark, p = pred_at_landmark),
          ggplot2::aes(
            x = .data[["t"]],
            y = .data[["p"]],
            colour = "Predicted value (landmark)"
          ),
          size = 4L,
          shape = 18L
        )
      }) +
      ggplot2::scale_y_continuous(
        name = dynamic_covariate,
        limits = c(y_min, y_max),
        sec.axis = ggplot2::sec_axis(
          ~ (. - y_min) / (y_max - y_min),
          name = "Survival probability",
          breaks = c(0, 0.25, 0.5, 0.75, 1),
          labels = c("0", "0.25", "0.50", "0.75", "1")
        )
      ) +
      ggplot2::scale_colour_manual(name = NULL, values = palette, breaks = names(palette)) +
      ggplot2::xlab("Time") +
      ggplot2::ggtitle(
        if (!is.null(cluster)) {
          paste0("Patient ", id, " (most likely cluster: ", cluster, ")")
        } else {
          paste("Patient", id)
        },
        subtitle = if (!is.null(probs)) {
          paste0(
            "Cluster membership probabilities: ",
            paste0(
              "cluster ",
              names(probs),
              " = ",
              sprintf("%.2f", probs),
              collapse = ", "
            )
          )
        }
      ) +
      ggplot2::theme_bw()

    myplot
  }
)
