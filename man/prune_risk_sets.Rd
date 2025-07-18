% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/landmarking.R
\name{prune_risk_sets}
\alias{prune_risk_sets}
\title{Prune a set of individuals from a risk set}
\usage{
prune_risk_sets(x, landmark, individuals)
}
\arguments{
\item{x}{An object of class \code{\link{LandmarkAnalysis}}.}

\item{landmark}{a landmark time}

\item{individuals}{Vector of individuals to be pruned from}
}
\value{
An object of class \code{\link{LandmarkAnalysis}} after having
pruned the individuals indicated in \code{individuals} from the risk set
at landmark time \code{landmark}.
}
\description{
Prune a set of individuals from a risk set
}
