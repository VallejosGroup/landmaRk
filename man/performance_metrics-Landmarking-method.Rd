% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/landmarking_performance.R
\name{performance_metrics,Landmarking-method}
\alias{performance_metrics,Landmarking-method}
\title{Performance metrics}
\usage{
\S4method{performance_metrics}{Landmarking}(x, landmarks, windows, c_index = TRUE, brier = TRUE)
}
\arguments{
\item{x}{An object of class \code{\link{Landmarking}}.}

\item{landmarks}{A numeric vector of landmark times.}

\item{windows}{A vector of prediction windows determining horizon times.}

\item{c_index}{A logical. If TRUE (default), C index is reported.}

\item{brier}{A logical. If TRUE (default), Brier score is reported.}
}
\description{
Computes concordance index (c-index) and Brier scores at the specified landmark
times and prediction windows.
}
