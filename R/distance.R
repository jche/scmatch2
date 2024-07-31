
# functions for working with the distance metric
library(tidyverse)



#' Coerce non-numeric factor columns to integer values
#'
#' @param covs 2-d dataframe of N x p, where N is number of data, p is
#'   dimension of features
#'
#' @return Coerced dataframe where factors are turned into integers.
#'   If these are scaled by K >> 1, and caliper is less than 1, then
#'   this will enforce exact matching on the category.
#'
#' @export
#'
#' @examples
#' df = data.frame( A = c( "A", "B", "A", "C" ) )
#' coerce_covs( df )
coerce_covs <- function(covs) {
  # TODO: Change this to using model.matrix to make dummy variables

  covs_coerced <-  covs %>%
    mutate(across(!where(is.numeric), ~as.numeric(as.factor(.))))
  return(covs_coerced)
}


#' Scale covariates by scaling factor specified by the "scaling" variable
#'
#' Note this converts a dataframe to a matrix.
#'
#' @param covs A dataframe of covariates
#' @param scaling A vector of length = ncol(covs), or a scalar
#'
#' @return Scaled matrix of same dimension as covs
#'
#' @export
#'
scale_covs <- function(covs, scaling){
  p <- ncol(covs)
  len_scaling <- length(scaling)

  if ((len_scaling != 1) & (len_scaling != p)){
    stop(paste0("Length of scaling must be 1 or ncol(covs). Your len_scaling is ",
                len_scaling, ", ncol(covs) is ", ncol(covs)))
  }

  if (len_scaling == 1) {
    scaling <- rep(scaling, p)
  }

  if (p > 1){
    covs_scaled <- covs %>%
      as.matrix() %*%
      diag(scaling)
  }else{
    covs_scaled <- covs %>%
      as.matrix() * scaling
  }

  return(covs_scaled)
}



#' Generate (#tx) by (#co) distance matrix
#'
#' @param df dataframe contains at least covariates and treatments
#' specified by "covs" and "treatment" argument
#' @param covs a vector of characters,
#' name of covariates to be matched; default all "X" variables
#' @param treatment name of the treatment variable; default "Z"
#' @param scaling vector of scaling constants for covariates
#' @param metric a character vector, distance metric to use.
#' Default is c("maximum", "euclidean", "manhattan")
#'
#' @return (#tx) by (#co) distance matrix
#'
#' @export
gen_dm <- function(df,
                   covs = get_x_vars(df),
                   treatment = "Z",
                   scaling = 1,
                   metric = c("maximum", "euclidean", "manhattan")) {

  # if (!all(covs %in% names(df))){
  #   missing_covs <- covs[!covs %in% names(df)]
  #   stop("covariates specified are not in the dataframe: ", paste(missing_covs,collapse=", "),call.=F )
  # }

  # if (!(treatment %in% names(df))){
  #   stop("Input treatment variable is wrong. You inputted: ", treatment )
  # }

  metric <- match.arg(metric)

  # pull out df with only covariates
  covs <- df %>%
    select(all_of(covs))

  covs_coerced <-  coerce_covs(covs)

  covs <- scale_covs(covs_coerced, scaling)

  # Subset tx/co covs correctly
  tx_obs <- which(as.logical(pull(df, {{treatment}})))
  co_obs <- which(as.logical(!pull(df, {{treatment}})))

  tx_covs <- matrix(covs[tx_obs,], ncol = ncol(covs))   # in case only one tx unit
  co_covs <- covs[co_obs,]

  # compute (#tx) x (#co) distance matrix
  dm <- flexclust::dist2(tx_covs,
                         co_covs,
                         method = metric)
  # Assign tx/co indices to row/col names of the dm
  rownames(dm) <- tx_obs
  colnames(dm) <- co_obs

  return(dm)
}






