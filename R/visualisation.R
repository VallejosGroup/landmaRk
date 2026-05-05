#' Plot longitudinal trajectory and predicted survival curve for one individual
#'
#' Produces a single-panel plot with a common time axis. The left of the
#' landmark dashed line shows the individual's observed longitudinal
#' measurements and model-fitted trajectory; the right shows their
#' predicted survival curve. The predicted value at the landmark that feeds
#' into the survival sub-model is highlighted. For LCMM fits, predicted
#' cluster probabilities are annotated on the plot.
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
#'
#' @examples
setMethod(
  "plot",
  "LandmarkAnalysis",
  function(x, id, landmark, dynamic_covariate, horizon = NULL, train = TRUE, ...) {

    # ---- Input validation ----
    if (missing(id) || is.null(id)) {
      stop("@id is required.")
    }
    if (missing(landmark) || !is.numeric(landmark) || length(landmark) != 1) {
      stop("@landmark must be a single numeric value.")
    }
    if (!(landmark %in% x@landmarks)) {
      stop(paste("Risk set for landmark time", landmark, "has not been computed."))
    }
    if (missing(dynamic_covariate) || is.null(dynamic_covariate)) {
      stop("@dynamic_covariate is required.")
    }
    predictions_slot <- if (train) x@survival_predictions else x@survival_predictions_test
    if (length(predictions_slot) == 0) {
      stop(if (train) {
        "No survival predictions found. Call predict_survival() first."
      } else {
        "No out-of-sample survival predictions found. Call predict_survival() with validation_fold > 0."
      })
    }

    # ---- Resolve horizon ----
    available <- grep(
      paste0("^", landmark, "-"),
      names(predictions_slot),
      value = TRUE
    )
    if (length(available) == 0) {
      stop(paste("No survival predictions available for landmark time", landmark))
    }
    if (is.null(horizon)) {
      if (length(available) == 1) {
        model_name <- available
        horizon <- as.numeric(sub(paste0(landmark, "-"), "", model_name))
      } else {
        stop(paste0(
          "Multiple horizons available for landmark ", landmark, ": ",
          paste(available, collapse = ", "),
          ". Please specify @horizon."
        ))
      }
    } else {
      model_name <- paste0(landmark, "-", horizon)
      if (!(model_name %in% names(predictions_slot))) {
        stop(paste(
          "No survival predictions for landmark", landmark, "and horizon", horizon
        ))
      }
    }

    # ---- Longitudinal observations ----
    obs_data <- x@data_dynamic[[dynamic_covariate]] |>
      dplyr::filter(.data[[x@ids]] == .env$id, .data[[x@times]] <= landmark)

    if (nrow(obs_data) == 0) {
      stop(paste(
        "No observations found for individual", id,
        "at or before landmark time", landmark
      ))
    }

    # ---- Predicted value at landmark (optional) ----
    long_preds <- if (train) {
      x@longitudinal_predictions[[as.character(landmark)]][[dynamic_covariate]]
    } else {
      x@longitudinal_predictions_test[[as.character(landmark)]][[dynamic_covariate]]
    }
    pred_at_landmark <- if (is.null(long_preds)) {
      NULL
    } else if (is.data.frame(long_preds)) {
      as.numeric(long_preds[as.character(id), 1L])
    } else {
      as.numeric(long_preds[as.character(id)])
    }

    # ---- Fitted trajectory on a dense time grid ----
    longitudinal_fit <- x@longitudinal_fits[[as.character(landmark)]][[
      dynamic_covariate
    ]]
    traj_df <- NULL
    if (!is.null(longitudinal_fit)) {
      t_start   <- min(obs_data[[x@times]], na.rm = TRUE)
      time_grid <- seq(t_start, landmark, length.out = 100L)
      newdata_grid <- stats::setNames(
        data.frame(id, time_grid),
        c(x@ids, x@times)
      ) |>
        dplyr::left_join(
          x@data_static |> dplyr::filter(.data[[x@ids]] == .env$id),
          by = x@ids
        )

      if (inherits(longitudinal_fit, "lmerMod")) {
        obs_data_full <- obs_data |>
          dplyr::left_join(
            x@data_static |> dplyr::filter(.data[[x@ids]] == .env$id),
            by = x@ids
          )
        pop_pred <- as.numeric(predict(
          longitudinal_fit,
          newdata = newdata_grid,
          re.form = NA
        ))
        ind_pred <- as.numeric(.predict_lme4(
          x            = longitudinal_fit,
          newdata      = newdata_grid,
          subject      = x@ids,
          test         = !train,
          newdata_long = obs_data_full
        ))
        traj_df <- data.frame(
          time       = rep(time_grid, 2L),
          prediction = c(pop_pred, ind_pred),
          type       = rep(
            c("Population average", "Individual"),
            each = length(time_grid)
          )
        )
      } else if (inherits(longitudinal_fit, "hlme")) {
        predicted_class <- longitudinal_fit$pprob |>
          dplyr::filter(.data[[x@ids]] == .env$id) |>
          dplyr::pull(class)
        if (length(predicted_class) == 0L) predicted_class <- 1L
        traj_df <- data.frame(
          time       = time_grid,
          prediction = as.numeric(
            lcmm::predictY(
              longitudinal_fit,
              newdata  = newdata_grid,
              var.time = x@times
            )$pred[, predicted_class]
          ),
          type = "Individual"
        )
      }
    }

    # ---- Individual-specific survival curve ----
    dataset <- if (train) x@survival_datasets[[model_name]] else x@survival_datasets_test[[model_name]]
    sf      <- predictions_slot[[model_name]]
    id_idx  <- which(dataset[[x@ids]] == id)
    if (length(id_idx) == 0L) {
      stop(paste("Individual", id, "not found in survival dataset for", model_name))
    }

    n_pts          <- 200L
    surv_times_rel <- seq(0, horizon - landmark, length.out = n_pts)
    sfit           <- summary(sf, times = surv_times_rel)
    surv_mat       <- matrix(sfit$surv, nrow = n_pts, ncol = nrow(dataset))
    surv_df        <- data.frame(
      time = surv_times_rel + landmark,
      surv = surv_mat[, id_idx]
    )

    # ---- Y-axis range (longitudinal scale with a small buffer) ----
    y_vals <- c(
      obs_data[[x@measurements]],
      if (!is.null(traj_df)) traj_df$prediction,
      if (!is.null(pred_at_landmark)) pred_at_landmark
    )
    y_raw <- range(y_vals, na.rm = TRUE)
    y_buf <- max(0.05 * diff(y_raw), 0.5)
    y_min <- y_raw[1L] - y_buf
    y_max <- y_raw[2L] + y_buf

    # Scale survival [0, 1] -> [y_min, y_max] for plotting on primary axis
    surv_df$surv_scaled <- surv_df$surv * (y_max - y_min) + y_min

    # ---- Assemble plot ----
    myplot <- ggplot2::ggplot() +
      # Survival curve
      ggplot2::geom_line(
        data = surv_df,
        ggplot2::aes(x = .data[["time"]], y = .data[["surv_scaled"]]),
        colour = "steelblue", linewidth = 0.8
      ) +
      # Observed longitudinal measurements
      ggplot2::geom_point(
        data = obs_data,
        ggplot2::aes(x = .data[[x@times]], y = .data[[x@measurements]])
      ) +
      # Vertical dashed line separating longitudinal from survival
      ggplot2::geom_vline(
        xintercept = landmark, linetype = "dashed", colour = "grey40"
      ) +
      # Predicted value at landmark (diamond, highlighted; only if available)
      (if (!is.null(pred_at_landmark)) {
        ggplot2::geom_point(
          data = data.frame(t = landmark, p = pred_at_landmark),
          ggplot2::aes(x = .data[["t"]], y = .data[["p"]]),
          colour = "tomato", size = 4L, shape = 18L
        )
      }) +
      ggplot2::scale_y_continuous(
        name   = dynamic_covariate,
        limits = c(y_min, y_max),
        sec.axis = ggplot2::sec_axis(
          ~ (. - y_min) / (y_max - y_min),
          name   = "Survival probability",
          breaks = c(0, 0.25, 0.5, 0.75, 1),
          labels = c("0", "0.25", "0.50", "0.75", "1")
        )
      ) +
      ggplot2::xlab("Time") +
      ggplot2::theme_bw()

    # Fitted longitudinal trajectory (lme4 / lcmm)
    if (!is.null(traj_df)) {
      myplot <- myplot +
        ggplot2::geom_line(
          data = traj_df,
          ggplot2::aes(
            x        = .data[["time"]],
            y        = .data[["prediction"]],
            linetype = .data[["type"]]
          ),
          colour = "tomato", linewidth = 0.8
        ) +
        ggplot2::scale_linetype_manual(
          name   = NULL,
          values = c("Individual" = "solid", "Population average" = "dashed")
        )
    } else {
      # LOCF: show a step function extended to the landmark so the carried-
      # forward value connects naturally to the highlighted diamond
      last_obs_val <- obs_data |>
        dplyr::arrange(.data[[x@times]]) |>
        dplyr::slice_tail(n = 1L) |>
        dplyr::pull(x@measurements)
      step_data <- obs_data |>
        dplyr::arrange(.data[[x@times]]) |>
        dplyr::select(dplyr::all_of(c(x@times, x@measurements))) |>
        dplyr::bind_rows(
          stats::setNames(
            data.frame(landmark, last_obs_val),
            c(x@times, x@measurements)
          )
        )
      myplot <- myplot +
        ggplot2::geom_step(
          data = step_data,
          ggplot2::aes(x = .data[[x@times]], y = .data[[x@measurements]]),
          colour = "tomato", linewidth = 0.8
        )
    }

    # LCMM cluster probability annotation
    if (
      !is.null(longitudinal_fit) &&
        inherits(longitudinal_fit, "hlme") &&
        longitudinal_fit$ng > 1L
    ) {
      pprob_row <- longitudinal_fit$pprob |>
        dplyr::filter(.data[[x@ids]] == .env$id) |>
        dplyr::select(dplyr::starts_with("prob"))
      if (nrow(pprob_row) > 0L) {
        prob_label <- paste(
          vapply(seq_len(ncol(pprob_row)), function(k) {
            sprintf("Class %d: %.1f%%", k, pprob_row[[1L, k]] * 100)
          }, character(1L)),
          collapse = "\n"
        )
        myplot <- myplot +
          ggplot2::annotate(
            "text",
            x      = landmark + 0.05 * (horizon - landmark),
            y      = y_max,
            label  = prob_label,
            hjust  = 0, vjust = 1, size = 3, colour = "grey30"
          )
      }
    }

    myplot
  }
)
