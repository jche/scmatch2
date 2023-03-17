
# functions for conducting matching



# wrapper function --------------------------------------------------------



#' Main matching function
#'
#' @param df 
#' @param metric 
#' @param caliper 
#' @param cal_method adaptive caliper, fixed caliper, only 1nn caliper
#' @param est_method 
#' @param return 
#' @param dist_scaling 
#' @param ... 
#'
#' @return df with a bunch of attributes
get_cal_matches <- function(df,
                            metric = c("maximum", "euclidean", "manhattan"),
                            caliper = 1,
                            cal_method = c("adaptive", "fixed", "1nn"),
                            est_method = c("scm", "scm_extrap", "average"),
                            return = c("sc_units", "agg_co_units", "all"),
                            dist_scaling = df %>%
                              summarize(across(starts_with("X"), 
                                               function(x) {
                                                 if (is.numeric(x)) 1/sd(x)
                                                 else 1000
                                               })),
                            ...) {
  metric <- match.arg(metric)
  cal_method <- match.arg(cal_method)
  est_method <- match.arg(est_method)
  return <- match.arg(return)
  args <- list(...)
  
  ### use cal_method: generate matches
  
  # get caliper matches
  scmatches <- df %>%
    gen_matches(
      scaling = dist_scaling, 
      metric = metric,
      caliper = caliper,
      cal_method = cal_method,
      ...)
  
  ### use est_method: scm or average
  
  scweights <- est_weights(df,
                           matched_gps = scmatches$matches,
                           dist_scaling = dist_scaling,
                           est_method = est_method,
                           metric = metric)
  
  ### use return: sc units, aggregated weights per control, all weights
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
  
  # aggregate results
  res <- m.data
  
  # store information about scaling and adaptive calipers
  adacalipers_df <- tibble(
    id = df %>% filter(Z == 1) %>% pull(id),
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
  
  return(res)
}


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






# matching method ---------------------------------------------------------

# get caliper matches
gen_matches <- function(df, 
                        covs = starts_with("X"),
                        treatment = Z,
                        scaling,
                        metric,
                        caliper,
                        cal_method = c("adaptive", "1nn", "fixed"),
                        ...) {
  cal_method <- match.arg(cal_method)
  args <- list(...)
  
  # store helpful constants
  ntx <- df %>% pull({{treatment}}) %>% sum()   # number of tx units
  p   <- df %>% select({{covs}}) %>% ncol()     # number of matched covariates
  
  ### step 0: generate distance matrix
  
  dm <- gen_dm(df, 
               covs={{covs}}, 
               treatment={{treatment}}, 
               scaling=scaling, 
               metric=metric)
  
  
  ### step 1: get caliper for each treated unit
  
  # adaptive caliper: min(caliper, 1nn)
  if (cal_method == "adaptive") {
    # default: p+1 nearest neighbors
    if (is.null(args$knn)) {
      knn <- p+1
    } else {
      knn <- args$knn
    }
    
    min_dists <- 1:ntx
    for (i in 1:ntx) {
      temp <- dm[i,]
      temp_sorted <- sort(temp)
      
      # idea:
      #  - if 1nn farther than caliper: caliper becomes 1nn
      #  - if 1nn closer than caliper: caliper becomes min(caliper, knn)
      if (temp_sorted[1] > caliper) {       # 1nn is farther than caliper
        min_dists[i] <- temp_sorted[1]
      } else {                              # 1nn is closer than caliper
        min_dists[i] <- min(caliper, temp_sorted[knn])
      }
    }
  } else if (cal_method == "1nn") {
    min_dists <- apply(dm, 1, min)
  } else {
    min_dists <- rep(caliper, nrow(dm))
  }
  
  
  ### step 2: remove all control units farther than caliper
  
  # per tx unit (row), remove all co units (cols) farther than min_dists away
  for (i in 1:ntx) {
    temp <- dm[i,]
    temp[temp > min_dists[i]] <- NA
    dm[i,] <- temp
    
    # if min_dists is too big (i.e., no exact match on discrete covs),
    #  ensure that no matches are returned
    if (min_dists[i] >= 1000) {
      dm[i,] <- NA
    }
  }
  
  ### step 3: generate df of matched controls for each treated unit
  df_list <- map(1:ntx, function(x) {
    
    # record which co units are matched to each tx unit
    #  - name: row number in df
    #  - value: col number in dm
    matched_obs <- which(!is.na(dm[x,]))
    
    # for each matched co unit, record distance from tx unit
    distances <- dm[x, matched_obs]
    matched_rows <- df[names(matched_obs),] %>%
      mutate(dist = distances)
    
    # if no matches, drop treated unit
    if (nrow(matched_rows) == 0) {
      # warning(glue::glue("Dropped treated unit {x} from data:
      #                       unable to exact-match categorical covariates"))
      return(NULL)
    }
    
    df %>%
      filter(as.logical({{treatment}})) %>%
      dplyr::slice(x) %>%
      mutate(dist = 0) %>%
      rbind(matched_rows) %>%
      mutate(subclass = x)
  })

  # store number of co matches per tx unit
  num_matches <- (!is.na(dm)) %>% 
    rowSums()
  
  return(list(matches = df_list %>% discard(is.null),   # drop unmatched tx units
              adacalipers = min_dists,
              dm = dm
              # num_matches = num_matches
              ))
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
    match_cols <- df %>%
      select(starts_with("X") & where(is.numeric)) %>%
      names()
    
    scweights <- map(matched_gps, 
                     ~gen_sc_weights(.x, match_cols, dist_scaling, metric),
                     .progress="Producing SCM units...")   # add progress bar
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







