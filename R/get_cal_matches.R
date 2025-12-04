

# Scale each by standard deviation, and exact match on binary or
# categorical covariates.
default_scaling <- function( data, covs ) {
  data %>%
    summarize(across(all_of(covs),
                     function(x) {
                       if (is.numeric(x)) 1/sd(x)
                       else 1000
                     }))
}



make_treatment_table <- function( matches ) {

  caliper = params(matches)$caliper
  treatment = params(matches)$treatment

  rs <- purrr::map_dfr( matches$matches,
                        function( blob ) {
                          if ( is.null( blob[["weights"]] ) ) {
                            blob$weights <- NA
                          }
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
    id = names( matches$adacalipers ),
    adacal = matches$adacalipers )

  if ( nrow(rs) > 0 ) {
    rs <- full_join( rs, adacalipers_df, by = "id" ) %>%
      mutate( nc = ifelse( is.na(nc), 0, nc ),
              feasible = ifelse( !is.na(adacal) & adacal <= caliper, 1, 0 ),
              matched = ifelse( nc > 0, 1, 0 ),
              id = as.character( id ) )
  } else {
    rs = adacalipers_df
    rs$nc = 0
    rs$ess = 0
    rs$max_dist = NA
    rs$feasible = 0
    rs$matched = 0
    rs$subclass = rs$id
  }

  rs
}





#' Caliper Synthetic Matching
#'
#' This is the core method of the CSM package. Conduct (adaptive)
#' radius matching with optional synthetic step on the resulting sets
#' of controls. The function creates a matched dataset by finding
#' control units that are similar to treatment units based on
#' specified covariates and distance metrics, then optionally weights
#' these controls using synthetic control methods.
#'
#' @inheritParams gen_matches
#'
#' @param data The data frame of data to be matched. Must contain
#'   treatment indicator and covariates for matching.
#' @param covs Specification of covariates to use for matching. Can be
#'   variable names, list of column numbers, or NULL.  If NULL, will
#'   use covariate names found in the scaling parameter, or default to
#'   all columns starting with "X".
#'
#' @param treatment Name of the column in \code{data} containing the
#'   treatment indicator (with values 0/1 or TRUE/FALSE).
#'
#' @param form Formula of form treatment ~ cov1 + cov2 + ...
#'   specifying the treatment variable and covariates to use for
#'   matching.  If not null, will override covs and treatment.
#' @param warn A logical indicating whether to warn about dropped
#'   units (those that cannot be matched within the caliper).
#' @param k Integer specifying the number of neighbors to use when
#'   \code{rad_method = "knn"}.
#'
#' @return An S3 object of class "csm_matches" containing:
#'   \itemize{
#'     \item \code{matches}: A list of data frames, each containing matched control units for a treated unit
#'     \item \code{adacalipers}: Vector of adaptive calipers for each treated unit
#'     \item \code{dm_trimmed}: Distance matrix with control units farther than caliper set to NA
#'     \item \code{dm_uncapped}: Original distance matrix without censoring
#'     \item \code{treatment_table}: Table of treated units with matching information
#'   }
#'   The object also has attributes storing the settings used for
#'   matching.
#'
#' @seealso \code{\link{gen_matches}} for the underlying matching
#'   function, \code{\link{est_weights}} for weight calculation,
#'   \code{\link{result_table}} for extracting results
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
#' # Convert to data frame
#' as.data.frame(mtch)
#'
#' @export
get_cal_matches <- function( data,
                             form = NULL,
                             covs = NULL,
                             treatment = NULL,
                             metric = c("maximum", "euclidean", "manhattan"),
                             caliper = 1,
                             rad_method = c("adaptive", "fixed", "1nn", "knn"),
                             est_method = c("scm", "scm_extrap", "average"),
                             scaling = NULL,
                             id_name = "id",
                             warn = TRUE ,
                             k = 1,
                             dm = NULL
) {
  metric <- match.arg(metric)
  rad_method <- match.arg(rad_method)
  est_method <- match.arg(est_method)

  # formula is of form treatment ~ covariates.
  # Extract these into 'treatment' and 'covs'
  if ( !is.null(form) ) {
    mf <- model.frame( form, data = data )
    treatment <- colnames( mf )[1]
    covs <- colnames( mf )[ -1 ]
  } else if ( is.null( covs ) && !is.null(scaling) && !is.null(names(scaling)) ) {
    # Pull covariates from scaling list
    covs <- names( scaling )
  } else if ( is.null( covs ) ) {
    # Default to all columns starting with "X"
    covs <- get_x_vars( data )
  }

  if ( is.null( treatment ) ) {
    warning( "treatment variable not specified; defaulting to 'Z'" )
    treatment <- "Z"
  }


  # Add ID column if not present
  if ( !( id_name %in% colnames(data) ) ) {
    n_zeros <- nchar(as.character(nrow(data)))
    data$id <- sprintf(paste0("U%0", n_zeros, "d"), 1:nrow(data))
  } else {
    data$id = as.character(data[[id_name]])
  }

  if ( is.null(scaling) ) {
    scaling <- default_scaling( data, covs )
  }

  scnames = names(scaling)
  if ( !is.null(scnames) ) {
    stopifnot( all( scnames == covs ) )
  }

  ### use rad_method: generate matches
  # get caliper matches
  scmatches <- gen_matches( data,
                            covs = covs,
                            treatment = treatment,
                            scaling = scaling,
                            metric = metric,
                            caliper = caliper,
                            rad_method = rad_method,
                            k = k,
                            id_name="id" )

  # scmatches$matches is a list of length ntx. Each element is a data
  # frame of matched controls for the given unit.


  ### use est_method to add individual unit weights: scm or average
  scweights <- est_weights( matched_gps = scmatches,
                            est_method = est_method)


  # Make a tibble of treated units with calipers and number of controls, etc.
  treatment_table <- make_treatment_table( scweights )

  # Possibly warn if treated units were lost
  unmatched_units <- filter( treatment_table, matched==0 )
  if ( warn && nrow(unmatched_units) > 0) {
    if ( nrow(unmatched_units) <= 20 ) {
      warning(glue::glue("Dropped the following treated units ({nrow(unmatched_units)} units from {nrow(treatment_table)}) from data: {paste(unmatched_units$id, collapse=\", \")}"))
    } else {
      warning(glue::glue("Dropped {nrow(unmatched_units)} of {nrow(treatment_table)} treated units from data.") )
    }
  }

  scmatches$matches <- scweights$matches
  scmatches$treatment_table <- treatment_table

  # keep a bunch of attributes around, just in case
  settings <- attr(scmatches, "settings" )

  settings$est_method = est_method

  attr( scmatches, "settings" ) <- settings

  class(scmatches) = "csm_matches"

  return(scmatches)
}




