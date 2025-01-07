
# The various methods for the csm_matches object



#' @return is.csm_matches: TRUE if object is a csm_matches object.
#'
#' @export
#'
#' @rdname csm_matches
#'
is.csm_matches <- function(x)
{
  inherits(x, "csm_matches")
}


#' @return `[`: pull out rows and columns of the dataframe.
#'
#' @rdname csm_matches
#' @export
`[.csm_matches` <- function(x, ...)
{
  as.data.frame(x)[...]
}




#' @return `[[`: pull out single element of dataframe.
#'
#' @rdname csm_matches
#' @export
`[[.csm_matches` <- function(x, ...)
{
  as.data.frame(x)[[...]]
}



#' @title Dimension of csm_matches object
#' @return dim: Dimension of csm_matches (as matrix)
#'
#' @rdname csm_matches
#' @export
#'
dim.csm_matches <- function(x, ... )
{
  return( dim( as.data.frame(x) ) )
}




#' @title Cast csm_matches to data.frame
#'
#' @param row.names NULL or a character vector giving the row names
#'   for the data frame.
#' @param optional logical. If TRUE, setting row names and converting
#'   column names is optional. (Currently ignored.)
#' @param ... additional arguments to be passed to the
#'   as.data.frame.list methods.  (Currently ignored.)
#'
#' @return as.data.frame: The matched data as a clean dataframe (no
#'   more attributes from csm_matches).
#' @rdname csm_matches
#'
#' @export
#'
as.data.frame.csm_matches <- function(
    x, row.names = NULL, optional = FALSE, return = "all", ...
) {
  rs <- switch(return,
               sc_units     = agg_sc_units(csm$matches),
               agg_co_units = agg_co_units(csm$matches),
               all          = bind_rows(csm$matches) )


  rs
}





#' Print method for csm_matches object
#' @param x object to print
#' @export
print.csm_matches <- function(x, ...) {
  ntx = length( x$matches )

  tt <- x$treatment_table
  if ( !is.null( tt ) ) {
    ntx_over = sum( tt$feasible != 1 )
    d_min = round( min( tt$max_dist ), digits = 3 )
    d_max = round( max( tt$max_dist ), digits = 3 )
  } else {
    ntx_over = sum( x$adacalipers > attr(x,"settings")$caliper, na.rm=TRUE )
    d_min = NA
    d_max = NA
  }
  nshow = min( ntx, 5 )

  settings <- attr(x,"settings")
  adas <- paste( round( x$adacalipers[1:nshow], 3 ), collapse = ", " )
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
    scaling = paste( round( scaling[1:5], digits=3 ), collapse = ", " )
    scaling = paste0( scaling, ", ..." )
  } else {
    scaling = paste( round( scaling, digits=3), collapse = ", " )
  }
  cal <- settings$cal
  cat( glue::glue( 'csm_matches: matching w/ {settings$metric} distance on {covs} \
                   {ntx} Treated units matched to control units ({ntx_over} above set caliper) \
                   Adaptive calipers: {adas} \
                   \tTarget caliper = {cal} \
                   \tMax distance ranges {d_min} - {d_max} \
                   scaling: {scaling}' ) )

  cat("\n")
}



#' Get full raw table of all the units
#'
#' Return dataframe with control units repeated, grouped with treated
#' unit with individual weights with those treated units.  Treated
#' units also included.
#'
#' @param csm A csm_matches object from a matching call.
#' @param feasible_only TRUE means only return treated units and
#'   matched controls for units that could be matched within the set
#'   caliper.
#' @param nonzero_weight_only TRUE means drop any control units with 0
#'   weight (e.g., due to scm weighting).
#' @return dataframe
#' @export
full_unit_table <- function( csm,
                             feasible_only = FALSE,
                             nonzero_weight_only = FALSE ) {
  rs <- csm$matches %>%
    bind_rows()

  if ( feasible_only ) {
    fs = csm$treatment_table %>%
      filter( feasible == 1 ) %>%
      pull( subclass )
    rs <- filter( rs, subclass %in% fs )
  }

  if ( nonzero_weight_only ) {
    rs <- filter( rs, weights > 0 )
  }

  return( rs )
}


#' Return table of calipers for all treated units
#'
#' @export
caliper_table <- function( csm ) {
  csm$treatment_table %>%
    dplyr::select( id, subclass, max_dist, adacal )
}


#' Return table of treated units
#' @export
feasible_units <- function( csm ) {
  csm$treatment_table %>%
    filter( feasible == 1 ) %>%
    dplyr::select( -feasible )
}

#' List all IDs of subclasses of units that are feasible
#'
#' For all units matched without adapting the caliper, get the
#' subclass IDs (i.e., the ids that link controls to treated units in
#' matched clusters).
#'
#' @return Numeric vector of IDs
#'
#' @export
feasible_unit_subclass <- function( csm ) {
  csm$treatment_table %>%
    filter( feasible == 1 ) %>%
    pull( subclass )
}


#' Obtain rows of treatment table for units that were not matched
#'
#' @return dataframe
#' @export
unmatched_units <- function( csm ) {
  csm$treatment_table %>%
    filter( matched == 0 ) %>%
    dplyr::select( -matched, -feasible )
}


#' Get table of aggregated results, with rows for synthetic units if
#' that was asked for
#'
#' @param return Possible values: NULL, "sc_units", "agg_co_units", or "all".
#'   How to aggregate units, if at all, in making the result table.
#' @param feasible_only TRUE means only return units which were given
#'   non-zero weight.
#'
#' @return dataframe
#' @export
#'
result_table <- function( csm,
                          return = c( "all", "sc_units", "agg_co_units" ),
                          feasible_only = FALSE ) {

  return = match.arg( return )

  # Swap result type if asked
  rs <- switch(return,
               sc_units     = agg_sc_units(csm$matches),
               agg_co_units = agg_co_units(csm$matches),
               all          = bind_rows(csm$matches) )



  if ( feasible_only ) {
    fs = csm$treatment_table %>%
      filter( feasible == 1 ) %>%
      pull( subclass )
    rs <- filter( rs, subclass %in% fs )
  }
  rs
}

