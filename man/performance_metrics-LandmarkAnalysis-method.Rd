% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance.R
\name{performance_metrics,LandmarkAnalysis-method}
\alias{performance_metrics,LandmarkAnalysis-method}
\title{Performance metrics}
\usage{
\S4method{performance_metrics}{LandmarkAnalysis}(x, landmarks, horizons, c_index = TRUE, brier = TRUE)
}
\arguments{
\item{x}{An object of class \code{\link{LandmarkAnalysis}}.}

\item{landmarks}{A numeric vector of landmark times.}

\item{horizons}{Vector of prediction horizons up to when the survival
submodel is fitted.}

\item{c_index}{A logical. If TRUE (default), C index is reported.}

\item{brier}{A logical. If TRUE (default), Brier score is reported.}
}
\description{
Computes concordance index (c-index) and Brier scores at the specified
landmark times and prediction horizons.
}
