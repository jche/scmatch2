kurtosis <- function(x){
  S_T = sd(x)
  kurt = mean( (x - mean(x))^4 ) / S_T^4
  return(kurt)
}

MCvar_SE <- function(x){
  S_T = sd(x); R <- length(x); k_T <- kurtosis(x)
  return(S_T^2 * sqrt( (k_T-1)/R ))
}

MCSE_SE <- function(x){
  return(sqrt(MCvar_SE(x)))
}

MCSE_bias <- function(x){
  S_T = sd(x); R <- length(x)
  return( sqrt(S_T^2 /R ))
}
