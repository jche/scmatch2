
# functions for conducting matching



# wrapper function --------------------------------------------------------

# TODO
get_cem_matches <- function(
    df,
    Z_FORMULA = as.formula(paste0("Z~",
                                  paste0(grep("^X", names(df), value=T),
                                         collapse="+"))),
    num_bins,
    est_method = c("average", "scm"),
    return = c("sc_units", "agg_co_units", "all")) {
  est_method <- match.arg(est_method)
  return <- match.arg(return)
  # print(Z_FORMULA)
  # print(head(df))
  m.out3 <- MatchIt::matchit(
    Z_FORMULA,
    data = df,
    method = "cem",
    estimand = "ATT",
    grouping = NULL,   # exact match factor covariates
    cutpoints = num_bins,
    k2k = F)   # not 1:1 matching
  m.data3 <- MatchIt::match.data(m.out3)
  # print(glue("Number of tx units kept: {sum(m.data3$Z)}"))

  # transform m.data3 into our format
  co_units <- m.data3 %>%
    filter(Z == 0)
  tx_units <- m.data3 %>%
    filter(Z == 1) %>%
    mutate(subclass_new = 1:n())  # give each tx unit its own subclass,
  # instead of letting multiple tx units share subclass

  # generate list of dfs, tx unit in first row
  # TODO: this is terribly slow, for no good reason...
  cem_matched_gps <- tx_units %>%
    split(f = ~subclass_new) %>%
    map(function(d) {
      d %>%
        bind_rows(
          co_units %>%
            filter(subclass == d$subclass) %>%
            mutate(subclass_new = d$subclass_new))
    }) %>%
    map(~.x %>%
          mutate(subclass = subclass_new) %>%
          select(-subclass_new, -weights))

  # estimate outcomes within cells
  scweights <- est_weights(df,
                           covs=names(df),
                           matched_gps = cem_matched_gps,
                           dist_scaling = df %>%
                             summarize(across(starts_with("X"),
                                              function(x) {
                                                if (is.numeric(x)) num_bins / (max(x) - min(x))
                                                else 1000
                                              })),
                           est_method = est_method,
                           metric = "maximum")

  # aggregate outcomes as desired
  m.data <- switch(return,
                   sc_units     = agg_sc_units(scweights),
                   agg_co_units = agg_co_units(scweights),
                   all          = bind_rows(scweights))

  unmatched_units <- setdiff(df %>% filter(Z==1) %>% pull(id),
                             m.data %>% filter(Z==1) %>% pull(id))
  if (length(unmatched_units) > 0) {
    warning(glue::glue("Dropped the following treated units from data:
                        \t {paste(unmatched_units, collapse=\", \")}"))
  }

  attr(m.data, "feasible_units") <- m.data %>%
    filter(Z) %>%
    pull(id)

  return(m.data)
}

# get the radius size for each treated unit
# input: dm, rad_method, ntx
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
get_matched_co_from_dm_trimmed <-
  function(df, dm_trimmed, treatment){
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

#
#'
#' Generate matches for each treatment over controls
#' using specified covariate names,
#' distance metric specified by the scaling parameter
#' and the method of metric.
#' Finally, there is a selection of the units
#'
#' Note: we allow controls to be repeatedly used
#'
#'
#' @param df A df containing at least covs and treatment
#' @param covs
#' @param treatment
#' @param scaling
#' @param metric Character specifying the distance metric.
#'  One of "maximum", "euclidean", "manhattan"
#' @param caliper Caliper on the scaled distance. Usually set to 1, and control matches
#' using the scaling paramter
#' @param rad_method
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
gen_matches <- function(df,
                        covs = starts_with("X"),
                        treatment = "Z",
                        scaling=1,
                        metric="maximum",
                        caliper=1,
                        rad_method = c("adaptive", "1nn", "fixed"),
                        ...) {
  ### Step -1: set up some constants
  rad_method <- match.arg(rad_method)
  args <- list(...)

  # store helpful constants
  ntx <- df %>% pull({{treatment}}) %>% sum()   # number of tx units
  p   <- df %>% select({{covs}}) %>% ncol()     # number of matched covariates

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
    get_matched_co_from_dm_trimmed(df, dm_trimmed,treatment)

  return(list(matches = df_list %>% discard(is.null),   # drop unmatched tx units
              adacalipers = radius_sizes,
              dm_trimmed = dm_trimmed,
              dm_uncapped = dm_uncapped))
}




# estimate within matched sets --------------------------------------------

#' Generate weights for list of matched sets
#'
#' @param df full dataframe
#' @param matched_gps list of matched sets: tx unit in first row
#' @param dist_scaling tibble with scaling for each covariate
#' @param est_method "scm" or "average"
#'
#' @return list of matched sets, with 'weights'
est_weights <- function(df,
                        covs,
                        matched_gps,
                        dist_scaling,
                        est_method = c("scm", "average"),
                        metric = c("maximum", "euclidean", "manhattan")) {
  est_method = match.arg(est_method)
  metric = match.arg(metric)

  if (est_method == "scm") {
    # generate SCM matching formula
    #  - NOTE: non-numeric columns crashed augsynth for some reason,
    #          so I still don't mess with them here.
    match_cols <- covs

    scweights <- map(matched_gps,
                     ~gen_sc_weights(.x, match_cols,
                                     dist_scaling,
                                     metric),
                     .progress="Producing SCM units...")
    # needs to modify gen_sc_weights to get clear:
    #   a) what type of matched_gps is required
    #   b) whether match_cols can be ignored
  } else if (est_method == "average") {
    scweights <- map(matched_gps,
                     function(x) {
                       x %>%
                         group_by(Z) %>%
                         mutate(weights = 1/n()) %>%
                         ungroup()
                     })
  }
  return(scweights)
}







