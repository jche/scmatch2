

# Test the various functionalities of the csm_matches_object


# The various methods for the csm_matches object




test_that("basic object functionality", {

  # Generate example data
  set.seed(4044440)
  dat <- gen_one_toy(nt = 5)

  # Perform matching
  mtch <- get_cal_matches(dat,
                          treatment = "Z",
                          metric = "maximum",
                          scaling = c(1/0.2, 1/0.2),
                          caliper = 0.25,
                          rad_method = "adaptive",
                          est_method = "scm")

  # View matching results
  mtch

  expect_true( is.csm_matches( mtch ) )
  expect_true( !is.csm_matches( as.data.frame( mtch ) ) )


  a = mtch[ 1, ]
  expect_true( is.data.frame( a ) )

  expect_true( mtch[[ 1, 4 ]] )

  dim(mtch)
  expect_true( is.numeric( dim( mtch )[1] ) )
  expect_true( is.numeric( dim( mtch )[2] ) )

  dd = as.data.frame( mtch )
  expect_true( is.data.frame( dd ) )
  expect_equal( nrow( dd ), nrow( result_table( mtch ) ) )

  expect_output( ss <- print( mtch ) )

  expect_output( ss <- summary( mtch ) )


  tb <- caliper_table( mtch )
  expect_true( is.data.frame( tb ) )
  expect_true( all( c("id","subclass","max_dist","adacal") %in% colnames( tb ) ) )
  expect_true( all( tb$id %in% mtch$treatment_table$id ) )

  ff = feasible_units( mtch )
  ff
  expect_true( all( ff$max_dist <= 0.25 ) )

  dl = tb %>%
    dplyr::filter( !( id %in% feasible_units( mtch )$id ) )
  dl
  expect_true( all( dl$max_dist > 0.25 ) )


  fs = feasible_unit_subclass( mtch )
  expect_true( all( fs %in% ff$id ) )

  um <- unmatched_units( mtch )
  expect_true( is.data.frame( um ) )
  expect_equal( nrow( um ), 0 )

  expect_warning(
    umt <- update_matches( mtch, dat, caliper = 0.35,
                           rad_method = "fixed",
                           est_method = "average" )
  )
  umt
  expect_true( is.csm_matches( umt ) )

  rs = result_table( umt )
  rs
  expect_true( sum( rs$Z ) < nrow( mtch$treatment_table ) )


  rs = result_table( umt, return = "exact" )
  rs
  expect_equal( nrow( rs ), 0 )

  rs = result_table( umt, return= "sc_units" )
  rs
  expect_equal( nrow(rs), 8 )

  rs = result_table( umt, return= "agg_co_units" ) %>%
    dplyr::filter( Z == 0 )
  rs
  expect_true( all(rs$id %in% as.character( dat$id ) ) )


  # And if we have no matches left?

  expect_warning(
    mtch <- get_cal_matches(dat,
                            form = Z ~ X1 + X2,
                            metric = "maximum",
                            scaling = c(1/0.2, 1/0.2),
                            caliper = 0.05,
                            rad_method = "fixed",
                            est_method = "scm")
  )
  expect_output(
    summary( mtch ),
    "0 treated units matched to 0 of 4 control units"
  )

  tt = CSM:::unmatched_units(mtch)
  expect_equal( nrow(tt), 5 )


})


#' Obtain table of aggregated results
#'
#' This method will give the controls as synthetic controls
#' (sc_units), the individual controls with repeat rows if controls
#' were used for different treated units (all), or aggregated controls
#' with only one row per control unit (agg_co_units).
#'
#' @param return How to aggregate units, if at all, in making the
#'   result table.  Possible values: "sc_units" (the synthetic control
#'   units), "agg_co_units" (the unique control units, with total
#'   weight across all their matches), or "all" (control units will be
#'   repeated if matched multiply). "exact" returns only exact matches
#'   (up to machine precision on distance). Defaults to "all".
#' @param feasible_only TRUE means only return units which were
#'   matched without expanding the caliper.
#'
#' @return dataframe of the treatment and control units.  This
#'   dataframe can be analyzed as an as-if experimental dataset.
#' @export
#'
result_table <- function( csm,
                          return = c( "all", "sc_units", "agg_co_units", "exact" ),
                          feasible_only = FALSE,
                          nonzero_weight_only = FALSE ) {

  return = match.arg( return )

  # Swap result type if asked
  rs <- switch(return,
               sc_units     = agg_sc_units(csm$matches),
               agg_co_units = agg_co_units(csm$matches),
               all          = bind_rows(csm$matches),
               exact        = bind_rows(csm$matches) %>%
                 filter( abs(dist) < 10*.Machine$double.eps )
  )

  if ( return == "exact" ) {
    rs <- rs %>%
      group_by( subclass ) %>%
      mutate( nc = n() - 1 ) %>%
      ungroup() %>%
      dplyr::filter( nc > 0 ) %>%
      dplyr::select( -nc )
  }

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
result_table <- function( csm,
                             feasible_only = FALSE,
                             nonzero_weight_only = FALSE ) {
  result_table( csm, "all", feasible_only = feasible_only,
                nonzero_weight_only = nonzero_weight_only )
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
update_matches <- function( res, ... ) {
  stopifnot( inherits( res, "csm_matches" ) )

  args = attr( res, "settings" )
  args = modifyList( args, list( ... ) )

  covs = attr( res, "covariates" )
  new_res <- get_cal_matches( data=test_df,
                              covs = covs,
                              treatment = args$treatment,
                              metric = args$metric,
                              caliper = args$caliper,
                              rad_method = args$rad_method,
                              est_method = args$est_method,
                              scaling = args$scaling,
                              id_name = args$id_name,
                              k = args$k )
  return( new_res )
}




