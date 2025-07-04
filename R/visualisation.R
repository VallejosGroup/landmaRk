#' Plots longitudinal trajectories and survival curves for landmarking models.
#'
#' @param x An object of class \code{\link{LandmarkAnalysis}}.
#' @param type A character taking the value \code{'survival'} (survival curves)
#'   or \code{'longitudinal'} (model trajectories of dynamic covariates).
#' @param id The identifier for the unit (subject) whose data will be plotted.
#' @param landmark Numeric indicating a landmark time
#' @param window Numeric indicating a prediction window
#' @param dynamic_covariate A character indicating a dynamic covariate
#' @param avg A logical (by default, \code{FALSE}) indicating whether LCMM
#'   predictions are conditioned on the predicted cluster (\code{avg = FALSE})
#'   or averaged across clusters (\code{avg = TRUE}). Ignored if the
#'   longitudinal model is not an LCMM.
#' @param ... Additional arguments passed to [survminer::ggadjustedcurves()]
#'   for plotting survival curves.
#'
#' @returns
#' @export
#'
#' @examples
setMethod(
  "plot",
  "LandmarkAnalysis",
  function(
    x,
    type = "survival",
    id = NULL,
    landmark = NULL,
    window = NULL,
    dynamic_covariate = NULL,
    avg = FALSE,
    ...
  ) {
    if (type == "survival") {
      if (length(x@survival_fits) == 0) {
        stop("Survival submodel has not been fitted.\n")
      } else if (is.null(landmark)) {
        stop("Argument @landmark is required when @type = 'survival'")
      } else if (is.null(dynamic_covariate)) {
        stop("Argument @window is required when @type = 'survival'")
      }
      if (!is.null(id)) {
        # If id is specified, generate survival plot for the relevant subject
        # Create the survival curve
        name <- paste0(landmark, "-", window)
        myplot <- survminer::ggsurvplot(
          survival::survfit(x@survival_fits[[name]]),
          data = x@survival_datasets[[name]],
          ...
        )

        # Modify X axis so survival plot starts at landmark time
        myplot$plot <- myplot$plot +
          ggplot2::scale_x_continuous(
            breaks = seq(from = 0, to = window, length.out = 4),
            labels = landmark +
              round(seq(from = 0, to = window, length.out = 4), 2)
          )

        myplot
      } else {
        # If id is not specified, generate survival plot for the survival model
        # Create list of plots that this function will return.
        plots <- list()
        for (name in names(x@survival_fits)) {
          # Retrieve the landmark time and prediction window for this sub-model fit.
          name_split <- unlist(strsplit(name, "-"))
          landmark <- as.numeric(name_split[1])
          window <- as.numeric(name_split[2])
          # Plot survival curve.
          plots[[name]] <- survminer::ggadjustedcurves(
            x@survival_fits[[name]],
            data = x@survival_datasets[[name]],
            ...
          ) +
            ggplot2::ggtitle(paste(
              "Landmark = ",
              landmark,
              "Prediction window = ",
              window
            )) +
            ggplot2::scale_x_continuous(
              breaks = seq(from = 0, to = window, length.out = 4),
              labels = landmark +
                round(seq(from = 0, to = window, length.out = 4), 2)
            )
        }
        plots
      }
    } else if (type == "longitudinal") {
      # id = NULL, landmark = NULL, window = NULL, dynamic_covariate = NULL, avg = FALSE
      if (is.null(id)) {
        stop("Argument @id is required when @type = 'longitudinal'.")
      } else if (is.null(landmark)) {
        stop("Argument @landmark is required when @type = 'longitudinal'")
      } else if (is.null(dynamic_covariate)) {
        stop(
          "Argument @dynamic_covariate is required when @type = 'longitudinal'"
        )
      }
      # Dataframe in long format with individuals at risk
      newdata <- .construct_data(
        x,
        dynamic_covariate = dynamic_covariate,
        x@risk_sets[[as.character(landmark)]],
        landmark
      )

      # Select the relevant subject
      newdata <- newdata |>
        filter(get(x@ids) == id)

      # Create plot with observations
      myplot <- newdata |>
        ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(get(x@times), get(x@measurements))) +
        ggplot2::xlab("Time") +
        ggplot2::ylab(dynamic_covariate)

      # Add model predictions
      longitudinal_fit <- x@longitudinal_fits[[as.character(landmark)]][[
        dynamic_covariate
      ]]
      if (inherits(longitudinal_fit, "lmerMod")) {
        # If model was fitted with lme4, use predict to make predictions
        predictions <- predict(longitudinal_fit, newdata = newdata)
        newdata[, "prediction"] <- predictions
        newdata[
          nrow(newdata) + 1,
          c("ParticipantId", "time", "prediction")
        ] <- c(
          id,
          landmark,
          x@longitudinal_predictions[[as.character(landmark)]][[
            dynamic_covariate
          ]][which(x@risk_sets[[as.character(landmark)]] == id)]
        )
        # Add predicted regression line to plot
        myplot <- myplot +
          ggplot2::geom_line(
            ggplot2::aes(get(x@times), !!sym("prediction")),
            data = newdata
          )
      } else if (inherits(longitudinal_fit, "hlme")) {
        # If model was fitted with lcmm, use predictY to make predictions
        predictions <- lcmm::predictY(longitudinal_fit, newdata = newdata)
        num_clusters <- ncol(predictions$pred)
        if (avg == TRUE) {
          # If averaging over cluster trajectories
          class_probabilities <- longitudinal_fit$pprob |>
            filter(get(x@ids) == id) |>
            select(starts_with("prob")) |>
            as.matrix()
          newdata[, "prediction"] <- predictions$pred %*% t(class_probabilities)
          newdata[
            nrow(newdata) + 1,
            c("ParticipantId", "time", "prediction")
          ] <- c(
            id,
            landmark,
            x@longitudinal_predictions[[as.character(landmark)]][[
              dynamic_covariate
            ]][as.character(id)]
          )
          # Add predicted regression line to plot
          myplot <- myplot +
            ggplot2::geom_line(
              ggplot2::aes(get(x@times), !!sym("prediction")),
              data = newdata
            )
        } else {
          # If predicting with the predicted cluster
          predicted_cluster <- longitudinal_fit$pprob |>
            filter(get(x@ids) == id) |>
            pull(class)
          newdata[, "prediction"] <- predictions$pred[, predicted_cluster]
          newdata[
            nrow(newdata) + 1,
            c("ParticipantId", "time", "prediction")
          ] <- c(
            id,
            landmark,
            x@longitudinal_predictions[[as.character(landmark)]][[
              dynamic_covariate
            ]][as.character(id)]
          )
          newdata[
            nrow(newdata) + 1,
            c("ParticipantId", "time", "prediction")
          ] <- c(
            id,
            landmark,
            x@longitudinal_predictions[[as.character(landmark)]][[
              dynamic_covariate
            ]][as.character(id)]
          )
          # Add predicted regression line to plot
          myplot <- myplot +
            ggplot2::geom_line(
              ggplot2::aes(get(x@times), !!sym("prediction")),
              data = newdata
            )
        }
      } else {
        stop("Only LCMM and LME are currently supported.")
      }

      myplot
    } else {
      stop(
        "Argument type must be survival or longitudinal.",
        "Additional options will be supported in the future.\n"
      )
    }
  }
)
