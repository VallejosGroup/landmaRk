% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival.R
\name{fit_survival}
\alias{fit_survival}
\title{Fits the specified survival model at the landmark times and up to the horizon
times specified by the user}
\usage{
fit_survival(
  x,
  formula,
  landmarks,
  horizons,
  method,
  dynamic_covariates = c(),
  include_clusters = FALSE
)
}
\arguments{
\item{x}{An object of class \code{\link{LandmarkAnalysis}}.}

\item{formula}{A formula to be used in survival sub-model fitting.}

\item{landmarks}{Numeric vector of landmark times.}

\item{horizons}{Vector of prediction horizons up to when the survival
submodel is fitted.}

\item{method}{Method for survival analysis, either "survfit" or "coxph".}

\item{dynamic_covariates}{Vector of time-varying covariates to be used
in the survival model.}

\item{include_clusters}{Boolean indicating whether to propagate cluster
membership to survival analysis.}
}
\value{
An object of class \code{\link{LandmarkAnalysis}}.
}
\description{
Fits the specified survival model at the landmark times and up to the horizon
times specified by the user
}
\details{
\subsection{Mathematical formulation}{

This function estimates the conditional probability of survival to horizon
\eqn{s+w}, conditioned on having survived to the landmark time, \eqn{s},
that is
\deqn{\pi_i(s+w \vert s) = P(T_i > s+w \vert T_i \ge s, \bar{x}_i(s)), }
where \eqn{i} denotes an individual's index, \eqn{T_i} is the time to event
outcome for individual \eqn{i} and \eqn{\bar{x}_i(s)} are the covariates
observed for individual \eqn{i}, including the observed history of dynamic
covariates.
}
}
