
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
  rs <- result_table( x, return=return )

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
  cat( glue::glue( 'csm_matches: matching with "{settings$metric}" distance\
                        match covariates: {covs} \
                   {ntx} Treated units matched to control units ({ntx_over} above set caliper) \
                   Adaptive calipers: {adas} \
                   \tTarget caliper = {cal} \
                   \tMax distance ranges {d_min} - {d_max} \
                   scaling: {scaling}' ) )

  cat("\n")
}


#' Summary method for csm_matches object
#' @param x object to summarize
#' @param ... Extra arguments (currently ignored).
#' @export
summary.csm_matches <- function(x, ...) {

  print.csm_matches(x)

  rs = result_table(x)
  rsC = filter( rs, Z == 0 )

  nunique1 = length( unique( rsC$id ) )
  nunique2 = length( unique( rsC$id[ rsC$weights > 0 ] ) )
  cat( glue::glue( "{nunique1} unique control units matched, {nunique2} with non-zero weight" ) )
  cat( "\n" )

  cat( "Treatment match pattern (before weighting):\n" )
  rs0 = filter( rs, Z == 0 )
  tb = as.numeric( table( rs$subclass ) )
  max = max( tb )
  tb[ tb > 7 ] = 7
  tb = table( tb )
  if ( "7" %in% names(tb) ) {
    w = which( names(tb) == "7" )
    names(tb)[w] = paste0( "7-", max )
  }
  print(tb)

  cat( "Treatment match pattern (after weighting):" )
  rs0 = filter( rs, weights > 0, Z == 0 )
  tb = table( table( rs0$subclass ) )
  print( tb )

  cat( "Control unit reuse:" )
  tb = table( table( rsC$id ) )
  print( tb )

  cat( "Summary of aggregated control weights\n" )
  rs2 = result_table( x, return = "agg_co_units" ) %>%
    filter( weights > 0 )
  print( summary( rs2$weights ) )
}



#' Return table of calipers for all treated units
#'
#' @export
caliper_table <- function( csm ) {
  csm$treatment_table %>%
    dplyr::select( id, subclass, max_dist, adacal )
}


#' Return table of treated units
#'
#' @return Dataframe of all the treated units (not controls) that were
#'   matched within a caliper.
#'
#' @export
feasible_units <- function( csm ) {
  csm$treatment_table %>%
    filter( feasible == 1 ) %>%
    dplyr::select( -feasible )
}



#' List all IDs of subclasses of units that are feasible
#'
#' For all units matched without expanding the caliper, get the
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


#' Obtain rows of treatment table for treated units that were not matched
#'
#' @return dataframe, one row per treated unit.
#'
#' @export
unmatched_units <- function( csm ) {
  csm$treatment_table %>%
    filter( matched == 0 ) %>%
    dplyr::select( -matched, -feasible )
}


#' Get table of aggregated results, with rows for synthetic units if
#' that was asked for
#'
#' This method will give the controls as synthetic controls
#' (sc_units), the individual controls with repeat rows if controls
#' were used for different treated units (all), or aggregated controls
#' with only one row per control unit (agg_co_units).
#'
#' @param return Possible values: "sc_units", "agg_co_units", or
#'   "all". How to aggregate units, if at all, in making the result
#'   table.  Defaults to "all".
#' @param feasible_only TRUE means only return units which were
#'   matched without expanding the caliper.
#'
#' @return dataframe
#' @export
#'
result_table <- function( csm,
                          return = c( "all", "sc_units", "agg_co_units" ),
                          feasible_only = FALSE,
                          nonzero_weight_only = FALSE ) {

  return = match.arg( return )

  # Swap result type if asked
  rs <- switch(return,
               sc_units     = agg_sc_units(csm$matches),
               agg_co_units = agg_co_units(csm$matches),
               all          = bind_rows(csm$matches) )



  if ( feasible_only ) {
    fs = feasible_unit_subclass(csm)
    rs <- filter( rs, subclass %in% fs )
  }

  if ( nonzero_weight_only ) {
    rs <- filter( rs, weights > 0 )
  }

  rs
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
#' @return dataframe, one row per treated unit and a row per control
#'   unit for each time it was used.
#' @seealso [result_table()]
#' @export
full_unit_table <- function( csm,
                             feasible_only = FALSE,
                             nonzero_weight_only = FALSE ) {
  result_table( csm, "all", feasible_only = feasible_only,
                nonzero_weight_only = nonzero_weight_only )
}



