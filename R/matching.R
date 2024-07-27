
# functions for conducting matching




#' get the radius size for each treated unit
#'
#' @param dm data
#' @param rad_method Method for radius
#' @param caliper Caliper size
#'
#' @export
get_radius_size <- function(dm,
                            rad_method,
                            caliper){
  ntx = nrow(dm)
  radius_sizes <- numeric(ntx)
  if (rad_method == "adaptive") {
    # Adaptive: r = max(caliper, 1nn)
    for (i in 1:ntx) {
      temp <- dm[i,]
      temp_sorted <- sort(temp)
      radius_sizes[i] <- max(caliper, temp_sorted[1])
    }
  } else if (rad_method == "1nn") {
    radius_sizes <- apply(dm, 1, min)
  } else {
    radius_sizes <- rep(caliper, nrow(dm))
  }
  return(radius_sizes)
}




# matching method ---------------------------------------------------------
# generate df of matched controls for each treated unit
get_matched_co_from_dm_trimmed <- function(df, dm_trimmed, treatment) {
  ntx <- nrow(dm_trimmed)
  df_trt <- df %>%
    filter(.data[[treatment]] == 1)
  for (x in 1:ntx){
    matched_obs <- which(!is.na(dm_trimmed[x,]))
  }

  map(1:ntx, function(x) {
    # Step 1: get which co units are matched to each tx unit
    matched_obs <- which(!is.na(dm_trimmed[x,]))
    # if no matches, drop treated unit
    if (length(matched_obs) == 0) {
      return(NULL)
    }


    # for each matched co unit, record distance from tx unit
    distances <- dm_trimmed[x, matched_obs]
    matched_rows <- df[names(matched_obs),] %>%
      ungroup() %>%
      mutate(dist = distances)

    ## Step 2: get the treat uni
    df_trt_x <- df_trt %>%
      ungroup() %>%
      slice(x)

    df_tmp <- df_trt_x %>%
      mutate(dist = 0) %>%
      rbind(matched_rows)

    df_tmp %>%
      mutate(subclass = x)
  })
}

set_NA_to_unmatched_co <- function(dm_uncapped, radius_sizes){
  ntx <- nrow(dm_uncapped)
  for (i in 1:ntx) {
    temp <- dm_uncapped[i,]
    temp[temp > radius_sizes[i]] <- NA
    dm_uncapped[i,] <- temp
  }
  return(dm_uncapped)
}






#' Match treated units to control units
#'
#' Generate matches for each treatment over controls using specified
#' covariate names, distance metric specified by the scaling parameter
#' and the method of metric. Finally, there is a selection of the
#' units
#'
#' Note: we allow controls to be repeatedly used (we match with
#' replacement).
#'
#'
#' @param df A data frame containing all the listed covariates and the
#'   treatment indicator.
#' @param covs List of covariates names.
#' @param treatment Variable in df that has treatment indicator (0/1
#'   variable).
#' @param scaling Scaling factor, defaults to 1.  Each covariate will
#'   be scaled by scaling factor.  Can be a vector of length P, with P
#'   being number of covariates.
#' @param metric Character specifying the distance metric. One of
#'   "maximum", "euclidean", "manhattan"
#' @param caliper Caliper on the scaled distance. Usually set to 1;
#'   control variable importance via the scaling paramter
#' @param rad_method Method to determine the adaptive radius size for
#'   each treated unit.  fixed means drop treated units with no
#'   matches closer than caliper. 1nn means largest of closest control
#'   and caliper.
#' @param ... Extra arguments
#'
#' @return A list of results.  matches: list of small datasets of
#'   matched control units for each treated unit.  adacalipers: vector
#'   of adaptive calipers for each treated unit.  dm_trimmed: distance
#'   matrix with control units farther than caliper censored with NA
#'   dm_uncapped: distance matrix without any control units censored.
#' @export
#'
#' @examples
#' df <- CSM:::gen_one_toy()
#' mtch <- gen_matches(df, covs = c("X1", "X2"), treatment = "Z")
#' names( mtch )
gen_matches <- function(df,
                        covs = get_x_vars(df),
                        treatment = "Z",
                        scaling=1,
                        metric="maximum",
                        caliper=1,
                        rad_method = c("adaptive", "1nn", "fixed"),
                        ...) {

  ### Step -1: set up some constants
  rad_method <- match.arg(rad_method)
  args <- list(...)

  stopifnot( treatment %in% colnames(df) )
  #stopifnot( all( covs %in% colnames(df) ) )

  # store helpful constants
  ntx <- df %>% pull({{treatment}}) %>% sum()   # number of tx units
  p   <- df %>% select(all_of(covs)) %>% ncol()     # number of matched covariates

  ### step 0: generate distance matrix, and
  #   store an uncapped version as well.
  dm_uncapped <- dm <- gen_dm(df,
                              covs=covs,
                              treatment=treatment,
                              scaling=scaling,
                              metric=metric)


  ### step 1: get the radius size for each treated unit
  radius_sizes <-
    get_radius_size(dm, rad_method, caliper)


  ### step 2: remove all control units farther than caliper
  # per tx unit (row), remove all co units (cols)
  # farther than radius_sizes away
  dm_trimmed <-
    set_NA_to_unmatched_co(dm_uncapped, radius_sizes)


  ### step 3: generate df of matched controls for
  #   each treated unit
  df_list <-
    get_matched_co_from_dm_trimmed(df, dm_trimmed, treatment)
  nulls = map_lgl(df_list, is.null)
  radius_sizes[nulls] = NA


  res <- list(matches = df_list %>% discard(is.null),   # drop unmatched tx units
              adacalipers = radius_sizes,
              dm_trimmed = dm_trimmed,
              dm_uncapped = dm_uncapped)

  class(res) <- "csm_matches"
  attr( res, "settings" ) <- list( caliper = caliper,
                                   rad_method = rad_method,
                                   metric = metric,
                                   scaling = scaling )
  attr( res, "covariates" ) <- covs
  return(res)
}



# calculate unit weights within matched sets --------------------------------------------

#' Estimate control weights within matched sets
#'
#' Generate weights for list of matched sets, typically by the
#' synthetic control method (but you can also simply average).
#'
#' @param matched_gps Either a csm_matches object or a list of matched
#'   sets, with a tx unit in first row of each set.
#' @param covs List of covariates to calculate similarity scores on if
#'   using scm.
#' @param scaling tibble with scaling for each covariate
#' @param est_method "scm" or "average"
#'
#' @return If passed an csm object, an updated csm object with
#'   weights. Otherwise, return list of matched sets, with 'weights'
#'   for each control unit.
#'
#' @export
est_weights <- function(matched_gps,
                        covs = attr( matched_gps, "covariates" ),
                        scaling = attr( matched_gps, "settings" )$scaling,
                        est_method = c("scm", "average"),
                        metric = attr( matched_gps, "settings" )$metric ) {
  #                          c("maximum", "euclidean", "manhattan")) {

  est_method = match.arg(est_method)
  metric = match.arg(metric, choices=c("maximum", "euclidean", "manhattan"))

  matches = matched_gps
  if ( is.csm_matches(matched_gps) ) {
    matches = matched_gps$matches
  }

  if (est_method == "scm") {
    # generate SCM matching formula
    #  - NOTE: non-numeric columns crashed augsynth for some reason,
    #          so I still don't mess with them here.
    match_cols <- covs

    scweights <- map( matches,
                      ~gen_sc_weights(.x, match_cols,
                                      scaling,
                                      metric),
                      .progress="Producing SCM units...")
    # needs to modify gen_sc_weights to get clear:
    #   a) what type of matched_gps is required
    #   b) whether match_cols can be ignored
  } else if (est_method == "average") {
    scweights <- map( matches,
                      function(x) {
                        x %>%
                          group_by(Z) %>%
                          mutate(weights = 1/n()) %>%
                          ungroup()
                      })
  }

  # We should not lose any units when calculating weights.
  stopifnot( length( matches ) == length( scweights ) )

  if ( is.csm_matches(matched_gps) ) {
    matched_gps$matches <- scweights
    return( matched_gps )
  } else {
    return(scweights)
  }
}







