

# Scale each by standard deviation, and exact match on binary
# covariates.
default_scaling <- function( df, covs ) {
  df %>%
    summarize(across(all_of(covs),
                     function(x) {
                       if (is.numeric(x)) 1/sd(x)
                       else 1000
                     }))
}



make_treatment_table <- function( df, matches, caliper ) {

  rs <- purrr::map_dfr( matches$matches,
                        function( blob ) {
                          summarise( blob,
                                     id = id[[1]],
                                     subclass = subclass[[1]],
                                     nc = n() - 1,
                                     ess = sum(weights[-1])^2 / sum( weights[-1]^2  ),
                                     max_dist = max(dist) )
                        }
  )

  # store information about scaling and adaptive calipers
  adacalipers_df <- tibble(
    id = df %>% filter(Z == 1) %>% pull(id),
    adacal = matches$adacalipers )

  rs <- left_join( rs, adacalipers_df, by = "id" ) %>%
    mutate( nc = ifelse( is.na(nc), 0, nc ),
            feasible = ifelse( !is.na(adacal) & adacal <= caliper, 1, 0 ),
            matched = ifelse( nc > 0, 1, 0 ),
            id = as.character( id ) )

  rs
}





#' Caliper Synthetic Matching
#'
#' Conduct (adaptive) radius matching with optional synthetic step on
#' the resulting sets of controls.
#'
#'
#' @param df The data frame to be matched.  Option id column to
#'   uniquely identify units.  If missing, will make a new ID column
#'   (with name `id`).
#' @param metric A string specifying the distance metric
#' @param caliper A numeric specifying the caliper size
#' @param rad_method adaptive caliper, fixed caliper, only 1nn caliper
#' @param est_method A string specifying the estimation method
#' @param return A string specifying what to return
#' @param scaling A vector of scaling constants for covariates (can
#'   also be a single row matrix).
#' @param warn A logical indicating whether to warn about dropped
#'   units.
#' @param id_name Name of column to look for ID values.  If that
#'   column not found, make canonical `id`.
#' @return df with a bunch of attributes.
#' @export
get_cal_matches <- function( df,
                             covs = starts_with("X"),
                             treatment = "Z",
                             metric = c("maximum", "euclidean", "manhattan"),
                             caliper = 1,
                             rad_method = c("adaptive", "fixed", "1nn","adaptive-5nn", "knn"),
                             est_method = c("scm", "scm_extrap", "average"),
                             return = c("sc_units", "agg_co_units", "all"),
                             scaling = default_scaling(df,covs),
                             id_name = "id",
                             warn = TRUE ,
                             k = 5  ## for KNN matching
                             ) {
  metric <- match.arg(metric)
  rad_method <- match.arg(rad_method)
  est_method <- match.arg(est_method)
  return <- match.arg(return)
  if ( !( id_name %in% colnames(df) ) ) {
    df$id <- paste0( "U", 1:nrow(df) )
  } else {
    df$id = as.character(df[[id_name]])
  }

  ### use rad_method: generate matches
  # get caliper matches
  scmatches <- df %>%
    gen_matches(
      covs = covs,
      treatment = treatment,
      scaling = scaling,
      metric = metric,
      caliper = caliper,
      rad_method = rad_method,
      k = k )

  # scmatches$matches is a list of length ntx.
  #   Each element is a data frame of matched controls


  ### use est_method to add individual unit weights: scm or average
  scweights <- est_weights( matched_gps = scmatches,
                            est_method = est_method)

  ### use return to figure out how to aggregate: sc units, aggregated weights per control, all weights
  res <- switch(return,
                sc_units     = agg_sc_units(scweights$matches),
                agg_co_units = agg_co_units(scweights$matches),
                all          = bind_rows(scweights$matches) )


  # Make a tibble of treated units with calipers and number of controls, etc.
  treatment_table <- make_treatment_table( df, scweights, caliper )

  # Possibly warn if treated units were lost
  unmatched_units <- filter( treatment_table, matched==0 )
  if ( warn && nrow(unmatched_units) > 0) {
    if ( length( nrow <= 20 ) ) {
      warning(glue::glue("Dropped the following treated units from data:
                        \t {paste(unmatched_units$id, collapse=\", \")}"))
    } else {
      warning(glue::glue("Dropped {nrow(unmatched_units) treated units from data.") )
    }
  }

  scmatches$matches <- scweights$matches
  scmatches$treatment_table <- treatment_table
  scmatches$result <- res

  # keep a bunch of attributes around, just in case
  settings <- attr(scmatches, "settings" )
  settings$rad_method = rad_method
  settings$est_method = est_method
  settings$return = return
  settings$treatment = treatment
  attr( scmatches, "settings" ) <- settings

  class(scmatches) = "csm_matches"

  return(scmatches)
}
