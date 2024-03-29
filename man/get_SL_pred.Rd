% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{get_SL_pred}
\alias{get_SL_pred}
\title{Get prediction from the fitted SuperLearner object}
\usage{
get_SL_pred(SL_fit, df_test, X_names)
}
\arguments{
\item{SL_fit}{A fitted SuperLearner object}

\item{X_names}{A vector of strings, specifying column names of covariates}

\item{df_pred}{The testing dataset}
}
\value{
A matrix of predictions. Column is the SL.lib. Rows are N_data
}
\description{
Get prediction from the fitted SuperLearner object
}
\examples{
mock_SL_fit <- create_mock_SL_fit()
mock_df_test <- create_mock_df_test()
X_names <- c("X1","X2")
result <- get_SL_pred(mock_SL_fit,
                      mock_df_test,
                      X_names)
}
