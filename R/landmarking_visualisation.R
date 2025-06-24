#' Plots survival curves for the fitted landmarking models.
#'
#' @param x An object of class Landmarking.
#' @param type A character string indicating the type of plot to generate.
#'   Currently, \code{type} must be "survival".  Additional options will be
#'   supported in the future.
#' @param ... Additional arguments passed to [survminer::ggadjustedcurves()]
#'   for plotting survival curves.
#'
#' @returns
#' @export
#'
#' @examples
setMethod(
  "plot",
  "Landmarking",
  function(x, type = "survival", ...) {
    if (type == "survival") {
      if (length(x@survival_fits) == 0) {
        stop("Survival submodel has not been fitted.\n")
      }
      # Create list of plots that this function will return.
      plots <- list()
      for (name in names(x@survival_fits)) {
        # Retrieve the landmark time and prediction window for this
        # sub-model fit.
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
            labels = landmark + seq(from = 0, to = window, length.out = 4)
          )
      }
      return(plots)
    } else {
      stop(
        "Argument type must be survival. Additional options will be ",
        "supported in the future.\n"
      )
    }
  }
)

plot_landmarking <- function(x, id, dynamic_covariate, landmark, window, ylim) {

  newdata <- x@data_dynamic[[dynamic_covariate]] |>
    filter(get(x@ids) == id) |>
    filter(get(x@times) <= landmark)

  landmark_prediction <-
    x@longitudinal_predictions[[as.character(landmark)]][[dynamic_covariate]][landmark]

  newdata$colour <- 1
  newdata <- rbind(newdata, c(id, landmark, landmark_prediction, 2))

  horizon <- landmark + window
  layout(
    mat = matrix(c(1,2), nrow = 1, ncol = 2),
    heights = c(1),
    widths = c(landmark/horizon, window/horizon)
  )
  # par(mfrow = c(1, 2))
  plot(
    newdata$time,
    newdata$value,
    col = newdata$colour,
    xlab = "Time",
    ylab = dynamic_covariate,
    pch = 16,
    ylim = ylim,
    xaxt = "n"
  )
  axis(1, at = round(seq(0, landmark, length.out = 5), 2))

  name <- paste0(landmark, "-", window)
  plot(
    survfit(x@survival_fits[[name]], data = x@data_static |> filter(get(x@ids) == id)),
    xaxt = "n",
    xlab = "Time",
    ylab = "Survival prob"
  )
  axis(1, at = seq(0, window, length.out = 4), labels = round(seq(landmark, landmark+window, length.out = 4), 2))
  #par(mfrow = c(1, 1))
}

plot_landmarking2 <- function(x, id, dynamic_covariate, landmark, window) {
  browser()
  cat("Hola")

  newdata <- x@data_dynamic[[dynamic_covariate]] |>
    filter(get(x@ids) == id) |>
    filter(get(x@times) <= landmark)

  landmark_prediction <-
    x@longitudinal_predictions[[as.character(landmark)]][[dynamic_covariate]][landmark]

  newdata$colour <- 1
  newdata <- rbind(newdata, c(id, landmark, landmark_prediction, 2))

  newdata |>
    ggplot(aes(x = get(x@times), y = get(x@measurements)), ) +
      geom_point(colour = newdata$colour)


  survminer::ggsurvplot(
    x@survival_fits[[name]],
    data = newdaa
  )
  name <- paste0(landmark, "-", window)
  survminer::ggadjustedcurves(
    survfit(x@survival_fits[[name]]),
    data = x@data_static |> filter(get(x@ids) == id)
  ) +
    ggplot2::ggtitle(paste(
      "Landmark = ",
      landmark,
      "Prediction window = ",
      window
    )) +
    ggplot2::scale_x_continuous(
      breaks = seq(from = 0, to = window, length.out = 4),
      labels = landmark + seq(from = 0, to = window, length.out = 4)
    )


  ggplot()
}
