
# functions for conducting matching




#' get the radius size for each treated unit
#'
#' @param dm Distance metric, rows are treated units, columns are control units
#' @param rad_method Method for determining radius size
#' @param caliper Caliper size
#' @param k For knn method, number of nearest neighbors
#'
#' @export
get_radius_size <- function(dm,
                            rad_method,
                            caliper,
                            k=1){
  ntx = nrow(dm)

  radius_sizes <- numeric(ntx)
  if (rad_method == "adaptive") {
    # Adaptive: r = max(caliper, knn)
    for (i in 1:ntx) {
      temp <- dm[i,]
      temp_sorted <- sort(temp)
      radius_sizes[i] <- max(caliper, temp_sorted[k])
    }
  } else if (rad_method == "1nn" || (rad_method=="knn" && k==1) ) {
    radius_sizes <- apply(dm, 1, min)
  } else if (rad_method == "knn") {
    for (i in 1:ntx) {
      temp <- dm[i,]
      temp_sorted <- sort(temp)
      radius_sizes[i] <- temp_sorted[k]
    }
  } else {
    radius_sizes <- rep(caliper, nrow(dm))
  }
  return(radius_sizes)
}




# matching method ---------------------------------------------------------

# generate dataframes of matched controls for each treated unit
get_matched_co_from_dm_trimmed <- function(data, dm_trimmed, treatment) {
  ntx <- nrow(dm_trimmed)
  data_trt <- data %>%
    filter(.data[[treatment]] == 1)

  # Add ID on the fly if needed, but try and use the one in the dataset.
  IDs <- paste0( "tx", 1:nrow(data_trt) )
  if ( "id" %in% colnames(data_trt) ) {
    IDs <- data_trt$id
  }

  # Unneeded I think
  #for (x in 1:ntx){
  #  matched_obs <- which(!is.na(dm_trimmed[x,]))
  #}

  map( 1:ntx, function(x) {
    # Step 1: get which co units are matched to each tx unit
    matched_obs <- which(!is.na(dm_trimmed[x,]))

    # if no matches, drop treated unit
    if (length(matched_obs) == 0) {
      return(NULL)
    }

    # for each matched co unit, record distance from tx unit
    distances <- dm_trimmed[x, matched_obs]
    matched_rows <- data[names(matched_obs),] %>%
      ungroup() %>%
      mutate(dist = distances)

    ## Step 2: get the treat unit
    data_trt_x <- data_trt %>%
      ungroup() %>%
      slice(x)

    data_tmp <- data_trt_x %>%
      mutate(dist = 0) %>%
      rbind(matched_rows)

    #   if ( nrow( data_tmp ) > 100 ) {
    #     browser()
    #   }

    data_tmp %>%
      mutate(subclass = IDs[[x]] )
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
#' Generally use the \code{get_cal_matches()} function instead, which
#' wraps this.
#'
#' @param metric Distance metric: "maximum", "euclidean", or
#'   "manhattan".
#'
#' @param caliper Caliper size. Sets the max allowed distance between
#'   matches. Often 1. Scaling usually handles covariate importance.
#'
#' @param rad_method How to set the radius for each treated unit:
#'   \itemize{
#'     \item "adaptive": max(caliper, distance to nearest control)
#'     \item "fixed": use caliper; drop treated with no controls inside it
#'     \item "1nn": distance to nearest neighbor
#'     \item "adaptive-5nn": adaptive with cap at 5th NN
#'     \item "knn": distance to the k-th nearest neighbor
#'   }
#'
#' @param est_method How to weight control units:
#'   \itemize{
#'     \item "scm": synthetic control weights
#'     \item "scm_extrap": SCM with extrapolation
#'     \item "average": simple average
#'   }
#'
#' @param scaling Length-P vector of scaling constants (or a
#'   single-row data frame). Each covariate is scaled by its value.
#'   Often 1/sd. Defaults to 1. If covs not supplied, the argument
#'   name is used to identify covariates.
#'
#' @param id_name Column name for unit IDs. If missing, an \code{id}
#'   column using row numbers is created.
#'
#' @param data Data frame with covariates and treatment indicator.
#'
#' @param covs Covariate names.
#'
#' @param treatment Treatment indicator column (0/1).
#'
#' @param k For adaptive, number of units needed for minimum radius.
#'   For "knn", number of neighbors. Defaults to 1.
#'
#' @param dm Optional treated-control distance matrix. Calipers are
#'   applied but distances are not recomputed.
#'
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
#' data <- CSM:::gen_one_toy()
#' mtch <- gen_matches(data, covs = c("X1", "X2"), treatment = "Z")
#' names( mtch )
gen_matches <- function( data,
                         covs = get_x_vars(data),
                         treatment = "Z",
                         scaling=1,
                         metric="maximum",
                         caliper=1,
                         rad_method = c("adaptive", "1nn", "fixed", "knn"),
                         id_name = NULL,
                         k = 1,
                         dm = NULL,
                         ...) {
  data <- data

  ### Step -1: set up some constants
  rad_method <- match.arg(rad_method)
  args <- list(...)

  if ( is.numeric(covs) ) {
    covs = colnames(data)[covs]
  }

  stopifnot( treatment %in% colnames(data) )
  #stopifnot( all( covs %in% colnames(data) ) )

  # store helpful constants
  ntx <- data %>% pull({{treatment}}) %>% sum()   # number of tx units
  p   <- data %>% dplyr::select(all_of(covs)) %>% ncol()     # number of matched covariates

  if ( !is.null( id_name ) ) {
    stopifnot( id_name %in% colnames(data) )
    data$id = data[[id_name]]
  } else {
    data <- data %>%
      group_by( across( all_of( treatment ) ) ) %>%
      mutate( id = paste0( ifelse( !.data[[treatment]], "co", "tx" ),
                           row_number() ) ) %>%
      ungroup()
  }
  data <- data %>%
    dplyr::relocate( id )

  ### step 0: generate distance matrix, and store an uncapped version
  #   as well (if it is not passed in)
  if ( is.null( dm ) ) {
    dm_uncapped <- dm <- gen_dm(data,
                                covs=covs,
                                treatment=treatment,
                                scaling=scaling,
                                metric=metric)
  } else {
    dm_uncapped <- dm
  }

  ### step 1: get the radius size for each treated unit
  radius_sizes <-
    get_radius_size(dm, rad_method, caliper, k)


  ### step 2: remove all control units farther than caliper
  #   per tx unit (row), remove all co units (cols)
  #.  farther than radius_sizes away
  dm_trimmed <-
    set_NA_to_unmatched_co(dm_uncapped, radius_sizes)


  ### step 3: generate dataframes of matched controls for
  #   each treated unit
  data_list <-
    get_matched_co_from_dm_trimmed(data, dm_trimmed, treatment)
  nulls = map_lgl(data_list, is.null)
  radius_sizes[nulls] = NA

  names( radius_sizes ) <- data %>%
    filter(.data[[treatment]] == 1) %>%
    pull(id)

  res <- list(matches = data_list %>% discard(is.null),   # drop unmatched tx units
              adacalipers = radius_sizes,
              dm_trimmed = dm_trimmed,
              dm_uncapped = dm_uncapped)

  class(res) <- "csm_matches"

  settings <- list( caliper = caliper,
                    rad_method = rad_method,
                    metric = metric,
                    scaling = scaling )
  settings$treatment = treatment
  settings$id_name = id_name
  settings$k = k

  attr( res, "settings" ) <- settings
  attr( res, "covariates" ) <- attr( dm, "covariates" )

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
est_weights <- function( matched_gps,
                         est_method = c("scm", "average"),
                         covs = params(matched_gps)$covariates,
                         treatment = params(matched_gps)$treatment,
                         scaling = params(matched_gps)$scaling,
                         metric = params(matched_gps)$metric ) {
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
                          group_by( across( all_of( treatment ) ) ) %>%
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







