#' S4 class for performing a landmarking analysis
#'
#' @slot landmarks A numeric vector of landmark times.
#' @slot data_static A data frame containing subject ids, static covariates,
#     censoring indicators  and event times.
#' @slot data_dynamic Data frame in long format with subject ids, measurement
#'  values, measurement times and measurement name.
#' @slot event_indicator Name of the column indicating event or censoring.
#' @slot ids Name of the column indicating subject ids.
#' @slot times Name of the column indicating observation time in
#'   \code{data_dynamic}.
#' @slot measurements Name of the column indicating measurement values in
#'   \code{data_dynamic}.
#' @slot event_time Name of the column indicating time of the event/censoring.
#' @slot risk_sets List of indices.
#' @slot longitudinal_fits List of model fits for the specified landmark times
#'   and biomarkers.
#' @slot longitudinal_predictions List of model predictions for the specified
#'   landmark times and biomarkers.
#' @slot survival_datasets List of survival dataframes used in the survival
#'   submodel.
#' @slot survival_fits List of survival model fits at each of the specified
#'   landmark times.
#' @slot survival_predictions List of time-to-event predictions for the
#'   specified landmark times and prediction horizons.
#'
#' @export
setClass(
  "LandmarkAnalysis",
  slots = c(
    landmarks = "numeric",
    data_static = "data.frame",
    data_dynamic = "list",
    event_indicator = "character",
    ids = "character",
    event_time = "character",
    times = "character",
    measurements = "character",
    risk_sets = "list",
    longitudinal_fits = "list",
    longitudinal_predictions = "list",
    survival_datasets = "list",
    survival_fits = "list",
    survival_predictions = "list"
  )
)

# Validator for object of class \code{\link{LandmarkAnalysis}}
#
# @param object An object of class \code{\link{LandmarkAnalysis}}.
#
# @returns TRUE if the input is valid, else a description of the problem
setValidity("LandmarkAnalysis", function(object) {
  error_str <- NULL
  if (is.null(names(object@data_dynamic))) {
    error_str <- c(
      error_str,
      "@data_dynamic must be a named list of dataframes"
    )
  }
  for (covariate in names(object@data_dynamic)) {
    if (!(object@ids %in% colnames(object@data_dynamic[[covariate]]))) {
      error_str <- c(
        error_str,
        "@ids must be a column in every dataframe in @data_dynamic"
      )
    }
    if (!(object@times %in% colnames(object@data_dynamic[[covariate]]))) {
      error_str <- c(
        error_str,
        "@times must be a column in every dataframe in @data_dynamic"
      )
    }
    if (
      !(object@measurements %in% colnames(object@data_dynamic[[covariate]]))
    ) {
      error_str <- c(
        error_str,
        "@measurements must be a column in every dataframe in @data_dynamic"
      )
    }
  }
  if (!(object@event_indicator %in% colnames(object@data_static))) {
    error_str <- c(
      error_str,
      "@event_indicator must be a column in dataframe @data_static"
    )
  }
  if (!(object@ids %in% colnames(object@data_static))) {
    error_str <- c(error_str, "@ids must be a column in dataframe @data_static")
  }
  if (!(object@event_time %in% colnames(object@data_static))) {
    error_str <- c(
      error_str,
      "@event_time must be a column in dataframe @data_static"
    )
  }
  if (length(error_str) == 0) {
    return(TRUE)
  } else {
    .eval_error_str(error_str)
  }
})

#' Creates an S4 class for a landmarking analysis
#'
#' @param data_static A data frame containing subject ids, static covariates,
#     censoring indicators  and event times.
#' @param data_dynamic Data frame in long format with subject ids, measurement
#'  values, measurement times and measurement name.
#' @param event_indicator  Name of the column indicating event or censoring.
#' @param ids Name of the column indicating subject ids.
#' @param event_time Name of the column indicating time of the event/censoring.
#' @param times Name of the column indicating observation time in
#'   \code{data_dynamic}.
#' @param measurements Name of the column indicating measurement values in
#'   \code{data_dynamic}.
#'
#' @returns An object of class \code{\link{LandmarkAnalysis}}
#' @export
LandmarkAnalysis <- function(
  data_static,
  data_dynamic,
  event_indicator,
  ids,
  event_time,
  times,
  measurements
) {
  # Find out static covariates of type characters
  char_columns <- names(which(sapply(data_static, class) == "character"))
  # Convert all character static covariates to factor, and raise a message
  if (length(char_columns) > 0) {
    data_static <- data_static |>
      mutate(across(where(is.character), as.factor))

    if (length(char_columns) == 1) {
      message(
        "Static covariate",
        paste(char_columns, collapse = ", "),
        " was coded as a character. Converted to factor."
      )
    } else {
      # More than one character column
      message(
        "Static covariates ",
        paste(char_columns, collapse = ", "),
        " were coded as characters. Converted to factors."
      )
    }
  }

  # Find out dynamic covariates of type characters
  for (dynamic_covariate in names(data_dynamic)) {
    if (measurements %in% names(data_dynamic[[dynamic_covariate]])) {
      if (
        inherits(data_dynamic[[dynamic_covariate]][, measurements], "character")
      ) {
        data_dynamic[[dynamic_covariate]][,
          measurements
        ] <- as.factor(data_dynamic[[dynamic_covariate]][, measurements])
        message(
          "Dynamic covariate ",
          dynamic_covariate,
          " were coded as character. Converted to factor."
        )
      }
    }
  }

  new(
    "LandmarkAnalysis",
    data_static = data_static,
    data_dynamic = data_dynamic,
    event_indicator = event_indicator,
    ids = ids,
    event_time = event_time,
    times = times,
    measurements = measurements,
    risk_sets = list(),
    longitudinal_fits = list(),
    longitudinal_predictions = list(),
    survival_fits = list(),
    survival_predictions = list()
  )
}

#' Displays an object of class "\code{\link{LandmarkAnalysis}}"
#'
#' @param object An object of class \code{\link{LandmarkAnalysis}}.
#'
#' @export
#'
#' @examples
setMethod(
  f = "show",
  signature = "LandmarkAnalysis",
  definition = function(object) {
    cat("Summary of LandmarkAnalysis Object:\n")
    cat("  Landmarks:", object@landmarks, "\n")
    cat("  Number of observations:", nrow(object@data_static), "\n")
    cat("  Event indicator:", object@event_indicator, "\n")
    cat("  Event time:", object@event_time, "\n")

    cat("  Risk sets:", "\n")
    if (length(object@risk_sets) > 0) {
      for (i in seq_along(object@risk_sets)) {
        n <- length(object@risk_sets[[i]])
        cat(
          "    Landmark ",
          object@landmarks[i],
          ": ",
          n,
          if (n == 1) " subject\n" else " subjects\n",
          sep = ""
        )
      }
    }
  }
)

#' Summarise a LandmarkAnalysis object
#'
#' @param object An object of class \code{\link{LandmarkAnalysis}}.
#' @param type If \code{longitudinal}, it summarises the longitudinal submodel.
#'  If \code{survival}, it summarises the survival submodel.
#' @param landmark A numeric indicating the landmark time.
#' @param horizon For survival submodels, a numeric indicating the horizon time.
#' @param dynamic_covariate For longitudinal submodels, a character indicating
#'  the dynamic covariate
#'
#' @returns A summary of the desired submodel
#' @export
#'
#' @examples
setMethod(
  "summary",
  signature(
    object = "LandmarkAnalysis"
  ),
  function(
    object,
    type = c("longitudinal", "survival"),
    landmark,
    horizon = NULL,
    dynamic_covariate = NULL
  ) {
    error_str <- NULL
    if (is.null(type) || !(type %in% c("longitudinal", "survival"))) {
      error_str <- c(error_str, "@type must be 'longitudinal' or 'survival'")
    }
    if (is.null(landmark) || !(is.numeric(landmark)) || length(landmark) > 1) {
      error_str <- c(error_str, "@landmark must be numeric")
    }

    .eval_error_str(error_str)

    if (type == "longitudinal") {
      # Summary of longitudinal submodel fit
      if (
        is.null(dynamic_covariate) ||
          !(is.character(dynamic_covariate)) ||
          length(dynamic_covariate) > 1
      ) {
        error_str <- c(
          error_str,
          "@dynamic_covariate must be of type character"
        )
      }
      if (!(as.character(landmark) %in% names(object@longitudinal_fits))) {
        stop(paste(
          "No longitudinal submodel has been fitted to landmark time",
          landmark
        ))
      } else if (
        !(dynamic_covariate %in%
          names(object@longitudinal_fits[[as.character(landmark)]]))
      ) {
        stop(paste(
          "No longitudinal submodel has been fitted for dynamic covariate",
          dynamic_covariate,
          "to landmark time",
          landmark
        ))
      }
      model_fit <- object@longitudinal_fits[[as.character(landmark)]][[
        dynamic_covariate
      ]]
      cat(
        capture.output(object@longitudinal_fits[[as.character(
          landmark
        )]][[dynamic_covariate]]),
        sep = "\n"
      )
    } else if (type == "survival") {
      # Summary of survival submodel fit
      if (is.null(horizon) || !(is.numeric(horizon)) || length(horizon) > 1) {
        error_str <- c(error_str, "@horizon must be numeric")
      }

      model_name <- paste0(landmark, "-", horizon)
      if (!(as.character(model_name) %in% names(object@survival_fits))) {
        error_str <- c(
          error_str,
          (paste0(
            "No survival submodel has been fitted to landmark-horizon time ",
            model_name,
            "\n"
          ))
        )
      }

      .eval_error_str(error_str)

      cat(capture.output(object@survival_fits[[model_name]]), sep = "\n")
    }
  }
)

# Accessor for landmarks
setGeneric("getLandmarks", function(object) standardGeneric("getLandmarks"))
setMethod("getLandmarks", "LandmarkAnalysis", function(object) object@landmarks)

# Accessor for event
setGeneric("getEvent", function(object) standardGeneric("getEvent"))
setMethod("getEvent", "LandmarkAnalysis", function(object) object@event)

# Accessor for ids
setGeneric("getIds", function(object) standardGeneric("getIds"))
setMethod("getIds", "LandmarkAnalysis", function(object) object@ids)

# Accessor for event_time
setGeneric("getEventTime", function(object) standardGeneric("getEventTime"))
setMethod("getEventTime", "LandmarkAnalysis", function(object) {
  object@event_time
})

# Accessor for risk_sets
setGeneric("getRiskSets", function(object) standardGeneric("getRiskSets"))
setMethod("getRiskSets", "LandmarkAnalysis", function(object) object@risk_sets)


#' Compute the list of individuals at risk at landmark times
#'
#' @param x An object of class \code{\link{LandmarkAnalysis}}.
#' @param landmarks Numeric vector of landmark times
#' @param .warn_when_less_than Integer indicating that a warning will be raised
#' when the number of observations prior to a landmark time is less than that
#' integer for certain individuals.
#' @param ... Additional arguments (not used)
#'
#' @returns An object of class \code{\link{LandmarkAnalysis}} including desired
#'   risk sets for the specified landmark times.
#' @export
#'
#' @details
#' A risk set describes all subjects still at risk (i.e., not experienced the
#' event of interest or censored) at a given time. In
#' \code{\link{LandmarkAnalysis}}, risk sets define which subjects should be
#' included in the longitudinal and survival sub-models for each landmark time.
#'
#' The risk sets are stored in the \code{risk_sets} slot of the
#' \code{\link{LandmarkAnalysis}} object, where each risk set is a list of
#' indices corresponding to the subjects at risk at the respective landmark
#' time.
#'
#' @examples
setGeneric(
  "compute_risk_sets",
  function(x, landmarks, .warn_when_less_than = 0, ...) {
    standardGeneric("compute_risk_sets")
  }
)


#' Compute the list of individuals at risk at landmark times
#'
#' @inheritParams compute_risk_sets
#'
#' @returns An object of class \code{\link{LandmarkAnalysis}}, including desired
#'   risk sets for the relevant landmark times.
#' @export
#'
#' @details
#' A risk set describes all subjects still at risk (i.e., not experienced the
#' event of interest or censored) at a given time. In \code{\link{LandmarkAnalysis}}, risk sets
#' define which subjects should be included in the longitudinal and survival
#' sub-models for each landmark time.
#'
#' The risk sets are stored in the \code{risk_sets} slot of the \code{\link{LandmarkAnalysis}}
#' object, where each risk set is a list of indices corresponding to the
#' subjects at risk at the respective landmark time.
#'
#' @examples
setMethod(
  "compute_risk_sets",
  "LandmarkAnalysis",
  function(x, landmarks, .warn_when_less_than = 0, ...) {
    `get(x@ids)` <- NULL # Global assignment to avoid R CMD check warning
    if (length(landmarks) == 1) {
      # If the vector of landmark times is of length 1
      if (landmarks %in% x@landmarks) {
        # Risk set for given landmark time is already in memory
        warning("Risk set for landmark time ", landmarks, " already computed")
      }
      # Add landmark time to the LandmarkAnalysis object
      x@landmarks <- c(x@landmarks, landmarks)
      # Compute risk set for given landmark time
      x@risk_sets[[as.character(landmarks)]] <-
        x@data_static[which(x@data_static[, x@event_time] >= landmarks), x@ids]

      # Now raise a warning if there are individuals with less than
      # @.warn_when_less_than observations prior to landmark time
      if (.warn_when_less_than > 1) {
        for (dynamic_covariate in names(x@data_dynamic)) {
          ids_few_observations <- x@data_dynamic[[dynamic_covariate]] |>
            # Filter observations prior to landmark time
            filter(.data[[x@times]] <= landmarks) |>
            # Group by individual id
            group_by(get(x@ids)) |>
            # Work out number of observations per individual
            summarise(n = n()) |>
            # Select individuals with less than @.warn_when_less_than observations
            filter(n < .warn_when_less_than) |>
            # Extract vector with individual ids
            pull(`get(x@ids)`)

          if (length(ids_few_observations) > 0) {
            warning(paste0(
              "The following individuals have less than ",
              .warn_when_less_than,
              " observations recorded prior to landmark time ",
              landmarks,
              " for dynamic covariate ",
              dynamic_covariate,
              ": ",
              paste(ids_few_observations, collapse = ", ")
            ))
          }
        }
      }
    } else {
      # Recursion to compute risk sets one-by-one
      x <- compute_risk_sets(x, landmarks[1], .warn_when_less_than)
      x <- compute_risk_sets(x, landmarks[-1], .warn_when_less_than)
    }
    x
  }
)

#' Prune a set of individuals from a risk set
#'
#' @param x An object of class \code{\link{LandmarkAnalysis}}.
#' @param landmark a landmark time
#' @param individuals Vector of individuals to be pruned from
#'
#' @returns An object of class \code{\link{LandmarkAnalysis}} after having
#'   pruned the individuals indicated in \code{individuals} from the risk set
#'   at landmark time \code{landmark}.
#' @export
#'
#' @examples
setGeneric(
  "prune_risk_sets",
  function(x, landmark, individuals) standardGeneric("prune_risk_sets")
)


#' Prune a set of individuals from a risk set
#'
#' @inheritParams prune_risk_sets
#'
#' @returns An object of class \code{\link{LandmarkAnalysis}} after having
#'   pruned the individuals indicated in \code{individuals} from the risk set
#'   at landmark time \code{landmark}.
#' @export
#'
#' @examples
setMethod(
  "prune_risk_sets",
  "LandmarkAnalysis",
  function(x, landmark, individuals) {
    error_str <- NULL

    if (!is.numeric(landmark) || length(landmark) != 1) {
      error_str <- c(
        error_str,
        "@landmark must be a numeric, indicating a landmark time"
      )
    } else if (!(as.character(landmark) %in% names(x@risk_sets))) {
      error_str <- c(
        error_str,
        paste("Risk set not available for landmark time", landmark)
      )
    }

    if (length(error_str) > 0) {
      .eval_error_str(error_str)
    }

    landmark <- as.character(landmark)

    # Calculate intersection between risk set and individuals to prune, and
    # work out the new risk set without those individuals.
    pruned_individuals <- intersect(x@risk_sets[[landmark]], individuals)

    if (length(pruned_individuals) > 0) {
      message(paste(
        "Removing",
        length(pruned_individuals),
        "individuals from risk set for landmark time",
        landmark
      ))
      x@risk_sets[[landmark]] <- setdiff(
        x@risk_sets[[landmark]],
        pruned_individuals
      )
    }

    # If any of the individuals to prune is not in risk, raise a warning
    if (length(pruned_individuals) < length(individuals)) {
      warning(paste(
        "A total of ",
        length(individuals) - length(pruned_individuals),
        "in @individuals are not in the risk set for landmark time",
        landmark
      ))
    }

    x
  }
)


# Accessor for survival fits
setGeneric(
  "getSurvivalFits",
  function(object) standardGeneric("getSurvivalFits")
)
setMethod(
  "getSurvivalFits",
  "LandmarkAnalysis",
  function(object) object@survival_fits
)

#' Prunes a landmark time from a \code{\link{LandmarkAnalysis}}, removing
#' the risk set, longitudinal submodel and survival submodel from the object.
#'
#' @param x An object of class \code{\link{LandmarkAnalysis}}.
#' @param ... Additional arguments
#'
#' @returns An object of class \code{\link{LandmarkAnalysis}}.
#' @export
#'
#' @examples
setGeneric(
  "prune",
  function(x, ...) standardGeneric("prune"),
  signature = "x"
)

#' Prunes a landmark time from a \code{\link{LandmarkAnalysis}}, removing
#' the risk set, longitudinal submodel and survival submodel from the object.
#'
#' @inheritParams prune
#' @param landmark A numeric indicating the landmark time.
#'
#' @returns An object of class \code{\link{LandmarkAnalysis}}.
#' @export
#'
#' @examples
setMethod(
  f = "prune",
  signature = "LandmarkAnalysis",
  definition = function(x, landmark = NULL) {
    if (is.null(landmark) || !is.numeric(landmark) || length(landmark) != 1) {
      stop("@landmark must be a numeric indicating a landmark time.")
    }
    landmark <- as.character(landmark)

    if (!(landmark %in% names(x@risk_sets))) {
      stop(paste0(
        "The risk set at landmark time ",
        landmark,
        " has not been computed."
      ))
    }
    message(paste0("Pruning landmark time ", landmark, "."))

    # Prune survival predictions
    model_names <- which(startsWith(
      names(x@survival_predictions),
      paste0(landmark, "-")
    ))
    if (length(model_names) > 0) {
      x@survival_predictions[[model_names]] <- NULL
    }

    # Prune survival model fits
    model_names <- which(startsWith(
      names(x@survival_fits),
      paste0(landmark, "-")
    ))
    if (model_names > 0) {
      x@survival_fits[[model_names]] <- NULL
    }

    # Prune survival datasets
    model_names <- which(startsWith(
      names(x@survival_datasets),
      paste0(landmark, "-")
    ))
    if (model_names > 0) {
      x@survival_datasets[[model_names]] <- NULL
    }

    # Prune longitudinal predictions
    model_name <- which(names(x@longitudinal_predictions) == landmark)
    if (model_name > 0) {
      x@longitudinal_predictions[[model_name]] <- NULL
    }

    # Prune longitudinal model fits
    model_name <- which(names(x@longitudinal_fits) == landmark)
    if (model_name > 0) {
      x@longitudinal_fits[[model_name]] <- NULL
    }

    # Prune risk sets
    x@risk_sets[[landmark]] <- NULL

    x
  }
)
