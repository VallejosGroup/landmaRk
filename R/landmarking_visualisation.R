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
