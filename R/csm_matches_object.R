


# The various methods for the csm_matches object




#' Check if object is a csm_matches object
#' @param x object to check
#' @return logical
#' @export
is.csm_matches <- function(x) {
  class(x) %in% "csm_matches"
}

#' Print method for csm_matches object
#' @param x object to print
#' @export
print.csm_matches <- function(x, ...) {
  ntx = length( x$matches )
  nshow = min( ntx, 5 )

  settings <- attr(x,"settings")
  adas <- paste( x$adacalipers[1:nshow], collapse = ", " )
  if ( nshow < ntx ) {
    adas = paste0( adas, ", ..." )
  }
  covs <- attr(x,"covariates")
  if ( length( covs ) > 5 ) {
    covs = paste( covs[1:5], collapse = ", " )
    covs = paste0( covs, ", ..." )
  } else {
    covs = paste( covs, collapse = ", " )
  }

  scaling = settings$scaling
  if ( length( scaling ) > 5 ) {
    scaling = paste( scaling[1:5], collapse = ", " )
    scaling = paste0( scaling, ", ..." )
  } else {
    scaling = paste( scaling, collapse = ", " )
  }
  cal <- settings$cal
  cat( glue::glue( 'csm_matches: matching w/ {settings$metric} distance on {covs} \
                   {ntx} Treated units matched to control units \
                   Adaptive calipers: {adas} \
                   \t(Target caliper = {cal}) \
                   scaling: {scaling}' ) )

  cat("\n")
}



#' @export
feasible_units <- function( csm ) {
  csm$treatment_table %>%
    filter( feasible == 1 ) %>%
    dplyr::select( -feasible )
}

#' @export
unmatched_units <- function( csm ) {
  csm$treatment_table %>%
    filter( matched == 0 ) %>%
    dplyr::select( -matched, -feasible )
}
