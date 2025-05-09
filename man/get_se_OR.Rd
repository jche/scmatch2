% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate.R
\name{get_se_OR}
\alias{get_se_OR}
\title{Get the standard error using the OR bootstrap approach}
\usage{
get_se_OR(
  matches,
  outcome = "Y",
  treatment = "Z",
  boot_mtd = "wild",
  B = 250,
  use_moving_block = F,
  seed_addition = 11,
  block_size = NA
)
}
\arguments{
\item{matches}{A CSM match object (R S3 object)}

\item{outcome}{Name of the outcome variable (default: "Y")}

\item{treatment}{Name of the treatment variable (default: "Z")}

\item{boot_mtd}{The bootstrap method to use. Options are "Bayesian", "wild", or "naive-resid".}

\item{B}{Number of bootstrap samples (default: 250)}

\item{use_moving_block}{Whether to use a moving block bootstrap (default: FALSE)}

\item{seed_addition}{Additional seed value to ensure reproducibility (default: 11)}

\item{block_size}{Block size for bootstrap (default: NA, automatically chosen)}
}
\value{
A tibble with standard error (SE), bootstrap confidence intervals, and sample sizes.
}
\description{
This function estimates the standard error of the ATT using a bootstrap approach
based on residuals from the OR method.
}

