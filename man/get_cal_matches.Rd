% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_cal_matches.R
\name{get_cal_matches}
\alias{get_cal_matches}
\title{Caliper Synthetic Matching}
\usage{
get_cal_matches(
  df,
  covs = starts_with("X"),
  treatment = "Z",
  metric = c("maximum", "euclidean", "manhattan"),
  caliper = 1,
  rad_method = c("adaptive", "fixed", "1nn", "adaptive-5nn"),
  est_method = c("scm", "scm_extrap", "average"),
  return = c("sc_units", "agg_co_units", "all"),
  scaling = default_scaling(df, covs),
  id_name = "id",
  warn = TRUE
)
}
\arguments{
\item{df}{The data frame to be matched.  Option id column to
uniquely identify units.  If missing, will make a new ID column
(with name \code{id}).}

\item{metric}{A string specifying the distance metric}

\item{caliper}{A numeric specifying the caliper size}

\item{rad_method}{adaptive caliper, fixed caliper, only 1nn caliper}

\item{est_method}{A string specifying the estimation method}

\item{return}{A string specifying what to return}

\item{scaling}{A vector of scaling constants for covariates (can
also be a single row matrix).}

\item{id_name}{Name of column to look for ID values.  If that
column not found, make canonical \code{id}.}

\item{warn}{A logical indicating whether to warn about dropped
units.}
}
\value{
df with a bunch of attributes.
}
\description{
Conduct (adaptive) radius matching with optional synthetic step on
the resulting sets of controls.
}
