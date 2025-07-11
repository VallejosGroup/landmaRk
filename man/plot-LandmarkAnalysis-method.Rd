% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualisation.R
\name{plot,LandmarkAnalysis-method}
\alias{plot,LandmarkAnalysis-method}
\title{Plots longitudinal trajectories and survival curves for landmarking models.}
\usage{
\S4method{plot}{LandmarkAnalysis}(
  x,
  type = "survival",
  id = NULL,
  landmark = NULL,
  window = NULL,
  dynamic_covariate = NULL,
  avg = FALSE,
  ...
)
}
\arguments{
\item{x}{An object of class \code{\link{LandmarkAnalysis}}.}

\item{type}{A character taking the value \code{'survival'} (survival curves)
or \code{'longitudinal'} (model trajectories of dynamic covariates).}

\item{id}{The identifier for the unit (subject) whose data will be plotted.}

\item{landmark}{Numeric indicating a landmark time}

\item{window}{Numeric indicating a prediction window}

\item{dynamic_covariate}{A character indicating a dynamic covariate}

\item{avg}{A logical (by default, \code{FALSE}) indicating whether LCMM
predictions are conditioned on the predicted cluster (\code{avg = FALSE})
or averaged across clusters (\code{avg = TRUE}). Ignored if the
longitudinal model is not an LCMM.}

\item{...}{Additional arguments passed to \code{\link[survminer:ggadjustedcurves]{survminer::ggadjustedcurves()}}
for plotting survival curves.}
}
\description{
Plots longitudinal trajectories and survival curves for landmarking models.
}
