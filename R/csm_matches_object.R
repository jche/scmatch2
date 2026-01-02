# R/csm_matches_object.R
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
  ntx = length( x$adacalipers )
  ntx_match = length( x$matches )

  tx_var = attr( x, "settings" )$treatment

  n_drop = sum( is.na( x$adacalipers ) )

  if ( length( x$matches ) > 0 ) {
    mtch = bind_rows( x$matches ) %>%
      dplyr::filter( .data[[tx_var]] == 0 )
    nco = n_distinct( mtch$id )

    mex <- mtch %>%
      filter( abs(dist) < 10*.Machine$double.eps )
    n_exact = n_distinct( mex$subclass )
  } else {
    nco = 0
    n_exact = 0
  }

  tt <- x$treatment_table
  if ( !is.null( tt ) ) {
    ntx_over = sum( tt$feasible != 1 & tt$matched == 1 )
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

  tco = dim( x$dm_uncapped )[[2]]

  if ( is.null( settings$est_method) ) {
    settings$est_method = "none"
  }
  cat( glue::glue( 'csm_matches: matching with "{settings$metric}" distance and "{settings$rad_method}" radii\
                        aggregating sets with "{settings$est_method}" method \
                        match covariates: {covs}' ) )
  cat("\n")
  cat( glue::glue( '{ntx_match} treated units matched to {nco} of {tco} control units \
                   \t({n_exact} exact matches, {ntx_match - ntx_over} below caliper, {ntx_over} above caliper) \
                   Adaptive calipers: {adas} \
                   \tTarget caliper = {cal} \
                   Max distance ranges {d_min} - {d_max} \
                   \tscaling: {scaling}' ) )

  cat("\n")
  if ( n_drop > 0 ) {
    cat( glue::glue( "{n_drop} treated units dropped" ) )
    cat("\n")
  }
}


#' Summary method for csm_matches object
#' @param x object to summarize
#' @param ... Extra arguments (currently ignored).
#' @export
summary.csm_matches <- function(x, outcome = NULL, ... ) {

  print.csm_matches(x)

  rs = result_table(x)

  if ( nrow(rs) == 0 ) {
    cat( "No matches were made.\n" )
    return( invisible(
      list( csm = x,
            att = NA,
            subclass_sizes = c( `0`= length(x$adacalipers) ) ) ) )
  }

  tx_var = attr( x, "settings" )$treatment
  rsC = filter( rs, .data[[tx_var]] == 0 )

  nunique1 = length( unique( rsC$id ) )
  nunique2 = length( unique( rsC$id[ rsC$weights > 0.000000001 ] ) )
  cat( glue::glue( "{nunique1} unique control units matched, {nunique2} with non-zero weight" ) )
  cat( "\n" )

  att = NULL
  if ( !is.null( outcome ) ) {
    stopifnot( outcome %in% names(x$matches[[1]]) )
    cat( "ATT estimates and sample sizes:\n" )
    att <- estimate_ATT(x, outcome=outcome) %>%
      dplyr::select( -V, -V_E, -V_P ) %>%
      as.data.frame()
    att %>%
      print( row.names = FALSE )
  }

  cat( "Subclass sizes (before weighting):\n" )
  tb = as.numeric( table( rsC$subclass ) )
  tbtb = table( tb )

  if ( length( tbtb ) > 7 ) {
    max = max( tb )

    # Get first existing value at least as big as 7 from tb
    first_val = as.numeric( names( tbtb )[[5]] )

    tb[ tb >= first_val ] = first_val

    tbtb = table( tb )
    w = which( names(tbtb) == as.character(first_val) )
    names(tbtb)[w] = paste0( first_val, "-", max )

  }

  print(tbtb)

  cat( "Subclass sizes (after weighting):" )
  rs0 = filter( rsC, weights > 0.000000001 )
  tb = table( table( rs0$subclass ) )
  print( tb )

  cat( "Control unit reuse (before weighting):" )
  tb = table( table( rsC$id ) )
  print( tb )

  cat( "Control unit reuse (after weighting):" )
  tb = table( table( rs0$id ) )
  print( tb )

  cat( "Summary of aggregated control weights\n" )
  rs2 = result_table( x, return = "agg_co_units" ) %>%
    filter( weights > 0.000000001 )
  print( summary( rs2$weights ) )

  return( invisible( list(
    csm = x,
    att = att,
    subclass_sizes = tbtb
  ) ) )
}



#' Return table of calipers for all treated units
#'
#' @export
caliper_table <- function( csm ) {
  csm$treatment_table %>%
    dplyr::select( id, subclass, feasible, min_dist, max_dist, adacal )
}


#' Return table of feasible treated units
#'
#' Return those treated units with controls within the set caliper.
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


#' Obtain table of treated units that were not matched
#'
#' @return dataframe, one row per treated unit.
#'
#' @export
unmatched_units <- function( csm ) {
  csm$treatment_table %>%
    filter( matched == 0 ) %>%
    dplyr::select( -matched, -feasible )
}


#' Obtain impact table
#'
#' This table summarizes the estimated impact for each treated unit,
#' along with the effective sample size of the controls used.
#'
#' @param scm A csm_matches object from a matching call.
#' @param outcome Name of the outcome variable.
#' @return dataframe with one row per treated unit, with columns:
#'   subclass, max_dist, outcome (estimated impact), precision
#'   (nominal precision of the impact estimate, calculated as 1/(1 +
#'   1/ess_C), nC (number of controls used).
#' @export
#'
#'
impact_table <- function( scm, outcome ) {

  tx_var = attr( scm, "settings" )$treatment

  rr = result_table( scm, include_caliper = TRUE, nonzero_weight_only = TRUE )

  rr$Y = ifelse( rr[[tx_var]], 1, -1 ) * rr[[ outcome ]]

  rsb <- rr %>%
    group_by( subclass, min_dist, max_dist ) %>%
    summarize( outcome = sum( Y * weights ),
               precision = 1 / (1 + 1/ess(weights[.data[[tx_var]]==0])),
               nC = n() - 1,
               .groups = "drop" )

  rsb
}



#' Obtain table of aggregated results
#'
#' This method will give the controls as synthetic controls
#' (sc_units), the individual controls with repeat rows if controls
#' were used for different treated units (all), or aggregated controls
#' with only one row per control unit (agg_co_units).
#'
#' @param csm A csm_matches object from a matching call.
#' @param feasible_only TRUE means only return treated units and
#'   matched controls for units that could be matched without
#'   expanding the caliper.
#' @param nonzero_weight_only TRUE means drop any control units with 0
#'   weight (e.g., due to scm weighting).
#' @param return How to aggregate units, if at all, in making the
#'   result table.  Possible values: "sc_units" (the synthetic control
#'   units), "agg_co_units" (the unique control units, with total
#'   weight across all their matches), or "all" (control units will be
#'   repeated if matched multiply). "exact" returns only exact matches
#'   (up to machine precision on distance). Defaults to "all".
#' @param include_caliper If TRUE, include columns for the caliper size
#'  and maximum distance for each treated unit.
#'
#' @return dataframe of the treatment and control units.  This
#'   dataframe can be analyzed as an as-if experimental dataset.
#'
#' @export
#'
result_table <- function( csm,
                          return = c( "all", "sc_units", "agg_co_units", "exact" ),
                          outcome = NULL,
                          id = NULL,
                          feasible_only = FALSE,
                          nonzero_weight_only = FALSE,
                          include_caliper = FALSE ) {

  return = match.arg( return )

  if ( !is.csm_matches( csm ) ) {
    stop( "Input csm must be a csm_matches object." )
  }

  keep = rep( TRUE, length( csm$matches ) )

  if ( feasible_only ) {
    keep[ csm$treatment_table$feasible != 1 ] <- FALSE
  }

  if ( !is.null( id ) ) {
    keep[ !csm$treatment_table$id %in% id ] <- FALSE
  }

  csm$matches <- csm$matches[ keep ]

  # Swap result type if asked
  rs <- switch(return,
               sc_units     = agg_sc_units(csm,outcome=outcome),
               agg_co_units = agg_co_units(csm,outcome=outcome),
               all          = bind_rows(csm$matches),
               exact        = bind_rows(csm$matches) %>%
                 filter( abs(dist) < 100*.Machine$double.eps )
  )

  # Deal with empty result table if nothing is matched.
  # (This is not great coding, I don't think.)
  if ( nrow(rs) == 0 ) {
    tx_var = attr( csm, "settings" )$treatment
    rs$id = character()
    rs$subclass = character()
    rs[tx_var] = integer()
    rs$weights = numeric()
  }


  if ( return == "exact" ) {
    # recalculate subgroup sizes and drop empty ones
    rs <- rs %>%
      group_by( subclass ) %>%
      mutate( nc = n() - 1 ) %>%
      ungroup() %>%
      dplyr::filter( nc > 0 ) %>%
      dplyr::select( -nc )
  }


  if ( nonzero_weight_only ) {
    rs <- filter( rs, weights > 0.000000001 )
  }

  if ( include_caliper ) {
    rs <- rs %>%
      left_join( caliper_table( csm ) %>%
                   dplyr::select( -id ),
                 by = "subclass" )
  }

  rs
}



# Backward compatible alias
full_unit_table <- function(csm,
                            feasible_only = FALSE,
                            nonzero_weight_only = FALSE,
                            ...) {
  result_table(
    csm,
    return = "all",
    feasible_only = feasible_only,
    nonzero_weight_only = nonzero_weight_only
  )
}



#' Update a matching call to change some parameters
#'
#' @param res A csm_matches object
#' @param ... Parameters to change in the matching call
#' @return A new csm_matches object with updated parameters
#'
#' @examples
#' # Generate example data
#' set.seed(4044440)
#' dat <- gen_one_toy(nt = 5)
#'
#' # Perform matching
#' mtch <- get_cal_matches(dat,
#'                         metric = "maximum",
#'                         scaling = c(1/0.2, 1/0.2),
#'                         caliper = 1,
#'                         rad_method = "adaptive",
#'                         est_method = "scm")
#'
#' # View matching results
#' mtch
#'
#' update_matches( mtch, caliper = 0.05, rad_method = "fixed" )
#'
#' @export
update_matches <- function( res, data, warn = TRUE, ... ) {
  stopifnot( inherits( res, "csm_matches" ) )

  args = attr( res, "settings" )
  args = modifyList( args, list( ... ) )

  covs = attr( res, "covariates" )
  new_res <- get_cal_matches( data=data,
                              covs = covs,
                              treatment = args$treatment,
                              metric = args$metric,
                              caliper = args$caliper,
                              rad_method = args$rad_method,
                              est_method = args$est_method,
                              scaling = args$scaling,
                              id_name = args$id_name,
                              warn = warn,
                              k = args$k )
  return( new_res )
}





#' Return the parameters of the method
#'
#' @param csm A csm_matches object
#' @return A list of the settings used in the matching call
#' @export
params <- function( csm ) {
  params = attr( csm, "settings" )
  params$covariates = attr( csm, "covariates" )

  params
}



#' Get the table of treated units
#'
#' @param csm A csm_matches object
#' @param id Optional list of ids for treated units
#' @param threshold Optional distance threshold for treated units;
#'   return only units below this threshold unless "bad" is set to
#'   TRUE.
#' @param bad If TRUE, return only treated units above the threshold.
#'   If bad is TRUE and threshold is NULL, return treated units above
#'   the set caliper.
#'
#' @return Dataframe of treated units
#'
#' @export
treatment_table <- function( csm, id = NULL, threshold = NULL,
                             bad = FALSE ) {
  tt <- csm$treatment_table
  if ( !is.null( id ) ) {
    tt <- tt %>%
      filter( id %in% id )
  }
  if ( bad || !is.null( threshold ) ) {
    if ( bad ) {
      if ( is.null( threshold ) ) {
        threshold = attr( csm, "settings" )$caliper
      }
      tt <- tt %>%
        filter( adacal > threshold )
    } else {
      tt <- tt %>%
        filter( adacal <= threshold )
    }
    }
  tt
}






#' Return result table for the bad matches
#'
#' @param csm A csm_matches object
#' @param threshold Distance threshold for bad matches
#' @param nonzero_weight_only TRUE means drop any control units with 0
#'   weight (e.g., due to scm weighting)
#'
#' @return Dataframe of bad matches
#' @export
bad_matches <- function( csm, threshold,
                         nonzero_weight_only = TRUE ) {
  rs <- result_table( csm, nonzero_weight_only = nonzero_weight_only  )

  tx_id = csm$treatment_table %>%
    filter( adacal > threshold ) %>%
    pull( id )
  tx_var = attr( csm, "settings" )$treatment

  bad_rs <- rs %>%
    filter( subclass %in% tx_id  ) %>%
    left_join( csm$treatment_table %>%
                 select( subclass, adacal ),
               by = "subclass"  ) %>%
    arrange( -adacal, subclass, -.data[[tx_var]] )

  bad_rs
}

