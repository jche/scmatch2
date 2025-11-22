

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
#' This is the core method of the CSM package. Conduct (adaptive)
#' radius matching with optional synthetic step on the resulting sets
#' of controls. The function creates a matched dataset by finding
#' control units that are similar to treatment units based on
#' specified covariates and distance metrics, then optionally weights
#' these controls using synthetic control methods.
#'
#' @param df The data frame of data to be matched. Must contain
#'   treatment indicator and covariates for matching.
#' @param covs Specification of covariates to use for matching. Can be
#'   variable names or tidyselect helpers like
#'   \code{starts_with("X")}. Defaults to variables starting with "X".
#'
#' @param treatment Name of the column in \code{df} containing the
#'   treatment indicator (0/1 or TRUE/FALSE). Defaults to "Z".
#'
#' @param form Formula of form treatment ~ cov1 + cov2 + ...
#'   specifying the treatment variable and covariates to use for
#'   matching.  If not null, will override covs and treatment.
#' @param metric A string specifying the distance metric. One of
#'   "maximum", "euclidean", or "manhattan".
#' @param caliper A numeric specifying the caliper size for matching.
#'   Controls the maximum allowable distance between matched units.
#' @param rad_method Method to determine the radius size for each
#'   treated unit. Options include:
#'   \itemize{
#'     \item "adaptive": Uses the maximum of the caliper and the distance to nearest control
#'     \item "fixed": Uses the specified caliper value
#'     \item "1nn": Uses the distance to the nearest neighbor
#'     \item "adaptive-5nn": Adaptive radius with cap at 5th nearest neighbor
#'     \item "knn": Uses the distance to the kth nearest neighbor
#'   }
#' @param est_method Method for estimating weights for control units.
#'   Options include:
#'   \itemize{
#'     \item "scm": Synthetic control method weighting
#'     \item "scm_extrap": Synthetic control with extrapolation
#'     \item "average": Simple average weighting
#'   }
#' @param scaling A vector of scaling constants for covariates (can
#'   also be a single row matrix). These are often just the inverses
#'   of covariate-specific standard deviations. Controls the relative
#'   importance of different covariates in the distance calculation.
#' @param id_name Name of column to look for ID values. If that column
#'   is not found, a canonical \code{id} corresponding to row numbers
#'   will be created.
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
#' function, \code{\link{est_weights}} for weight calculation,
#' \code{\link{result_table}} for extracting results
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
get_cal_matches <- function( df,
                             form = NULL,
                             covs = starts_with("X"),
                             treatment = "Z",
                             metric = c("maximum", "euclidean", "manhattan"),
                             caliper = 1,
                             rad_method = c("adaptive", "fixed", "1nn","adaptive-5nn", "knn"),
                             est_method = c("scm", "scm_extrap", "average"),
                             scaling = default_scaling(df,covs),
                             id_name = "id",
                             warn = TRUE ,
                             k = 5  ## for KNN matching
) {
  metric <- match.arg(metric)
  rad_method <- match.arg(rad_method)
  est_method <- match.arg(est_method)

  # formula is of form treatment ~ covariates.
  # Extract these into 'treatment' and 'covs'
  if ( !is.null(form) ) {
    mf <- model.frame( form, data = df )
    treatment <- colnames( mf )[1]
    covs <- colnames( mf )[ -1 ]
  }

  # Add ID column if not present
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
      k = k,
      id_name="id" )

  # scmatches$matches is a list of length ntx. Each element is a data
  # frame of matched controls for the given unit.


  ### use est_method to add individual unit weights: scm or average
  scweights <- est_weights( matched_gps = scmatches,
                            est_method = est_method)


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

  # keep a bunch of attributes around, just in case
  settings <- attr(scmatches, "settings" )

  settings$covariates = covs
  settings$treatment = treatment
  settings$metric = metric
  settings$caliper = caliper
  settings$rad_method = rad_method
  settings$est_method = est_method
  settings$scaling = scaling
  settings$id_name = id_name
  settings$k = k

  attr( scmatches, "settings" ) <- settings

  class(scmatches) = "csm_matches"

  return(scmatches)
}




