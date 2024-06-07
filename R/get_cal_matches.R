

#' Caliper Synthetic Matching
#'
#' This function implements (adaptive) radius matching with optional
#' synthetic step on the resulting sets of controls.
#'
#' @param df The data frame to be matched
#' @param metric A string specifying the distance metric
#' @param caliper A numeric specifying the caliper size
#' @param rad_method adaptive caliper, fixed caliper, only 1nn caliper
#' @param est_method A string specifying the estimation method
#' @param return A string specifying what to return
#' @param dist_scaling A vector of scaling constants for covariates
#' @param ... Additional arguments
#'
#' @return df with a bunch of attributes.
#' @export
get_cal_matches <- function( df,
                             covs = starts_with("X"),
                             treatment = "Z",
                             metric = c("maximum", "euclidean", "manhattan"),
                             caliper = 1,
                             rad_method = c("adaptive", "fixed", "1nn"),
                             est_method = c("scm", "scm_extrap", "average"),
                             return = c("sc_units", "agg_co_units", "all"),
                             dist_scaling = df %>%
                               summarize(across(all_of(covs),
                                                function(x) {
                                                  if (is.numeric(x)) 1/sd(x)
                                                  else 1000
                                                })),
                             ...) {
  metric <- match.arg(metric)
  rad_method <- match.arg(rad_method)
  est_method <- match.arg(est_method)
  return <- match.arg(return)
  args <- list(...)
  df$id <- 1:nrow(df)

  ### use rad_method: generate matches
  # get caliper matches
  scmatches <- df %>%
    gen_matches(
      covs = covs,
      treatment = treatment,
      scaling = dist_scaling,
      metric = metric,
      caliper = caliper,
      rad_method = rad_method,
      ...)

  # scmatches$matches is a list of length ntx.
  #   Each element is a data frame of matched controls

  ### use est_method: scm or average
  scweights <- est_weights(df,
                           covs=covs,
                           matched_gps = scmatches$matches,
                           dist_scaling = dist_scaling,
                           est_method = est_method,
                           metric = metric)

  ### use return: sc units, aggregated weights per control, all weights
  m.data <- switch(return,
                   sc_units     = agg_sc_units(scweights),
                   agg_co_units = agg_co_units(scweights),
                   all          = bind_rows(scweights))

  unmatched_units <- setdiff(df %>% filter(.data[[treatment]]==1) %>% pull(id),
                             m.data %>% filter(.data[[treatment]]==1) %>% pull(id))
  if (length(unmatched_units) > 0) {
    if ( length( unmatched_units <= 20 ) ) {
      warning(glue::glue("Dropped the following treated units from data:
                        \t {paste(unmatched_units, collapse=\", \")}"))
    } else {
      warning(glue::glue("Dropped {length(unmatched_units) treated units from data.") )
    }
  }

  # aggregate results
  res <- m.data

  # store information about scaling and adaptive calipers
  adacalipers_df <- tibble(
    id = df %>% filter(Z == 1) %>% pull(id) %>% as.character(),
    adacal = scmatches$adacalipers)
  attr(res, "scaling") <- dist_scaling
  attr(res, "adacalipers") <- adacalipers_df

  # store information about feasible units/subclasses
  feasible_units <- adacalipers_df %>%
    filter(adacal <= caliper) %>%
    pull(id)
  attr(res, "unmatched_units") <- unmatched_units
  attr(res, "feasible_units")  <- feasible_units
  attr(res, "feasible_subclasses") <- res %>%
    filter(id %in% feasible_units) %>%
    pull(subclass)

  # keep a bunch of data around, just in case
  attr(res, "scweights") <- scweights
  attr(res, "dm") <- scmatches$dm
  attr(res, "dm_uncapped") <- scmatches$dm_uncapped

  return(res)
}
