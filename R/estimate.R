


# functions for estimating effects
# The following two functions should be combined
#   once we get to know the input of get_att_ests
get_est_att_from_wt <- function(df,
                                input_wt){
  df_est_att <- df %>%
    cbind(wt=input_wt) %>%
    group_by(Z) %>%
    summarise(Y_wtd = weighted.mean(Y,wt))

  est_att <- diff(df_est_att$Y_wtd)
  return(est_att)
}


#' Estimate the ATT from a matched dataframe
#'
#' Given a matched dataset, calculate the estimated ATT
#'
#' @param matched_df A matched dataset
#'
#' @export
get_att_ests <- function(matched_df, outcome = "Y") {

  if ( is.csm_matches(matched_df) ) {
    matched_df <- matched_df$result
  }
  stopifnot( all( c("Z", "Y") %in% names(matched_df) ) )
  matched_df$Y = matched_df[[outcome]]

  matched_df %>%
    group_by(Z) %>%
    summarize(mn = sum(Y*weights) / sum(weights)) %>%
    summarize(est = last(mn) - first(mn)) %>%
    pull(est)
}

# get_att_ests_from_wt<-
#   function(df, input_wt)



calc_N_T_N_C <- function(preds_csm){
  if ( is.csm_matches( preds_csm ) ) {
    preds_csm <- full_unit_table(preds_csm)
  }

  N_T <- nrow(preds_csm %>% filter(Z==T))
  tmp <- preds_csm %>%
    filter(Z==F) %>%
    group_by(id) %>%
    summarise(w_i = sum(weights), .groups="drop")
  N_C_tilde <- N_T^2 / sum(tmp$w_i^2)
  return(list(N_T = N_T,
              N_C_tilde = N_C_tilde ))
}


ess <- function( weights ) {
  sum( weights )^2  / sum( weights^2 )
}

#' Calculate the variance for each subclass
#'
#' @param matches_filtered The filtered matched object
#' @param outcome Name of the outcome variable (default "Y")
#' @return A data frame with subclass variances
calculate_subclass_variances <-
  function(matches_filtered, outcome = "Y") {
  matches_filtered %>%
    group_by(subclass) %>%
    summarise(
      nj = n(),
      w_nj = ess(weights),
      var_cluster = var(!!sym(outcome)),
      .groups = "drop"
    )
}

#' Calculate the weighted average variance
#' NOTE: Could weight by nj, w_nj, or 1.  w_nj takes into account how
# much uncertainty is in each group, and thus might perform better
# with heteroskedasticity?
#'
#' @param cluster_var_df Data frame. Each row is a subclass and
#'  its variances, its weight
#' @param var_weight_type
#' @return Weighted average of subclass variances
calculate_weighted_variance <-
  function(cluster_var_df,
           var_weight_type = "num_units") {
  if (var_weight_type == "num_units"){ # number of units in the subclass
    weight = cluster_var_df$nj
  }else if (var_weight_type == "ess_units"){ # effective size of units in the subclass
    weight = cluster_var_df$w_nj
  }else if (var_weight_type == "uniform"){
    weight = rep(1, nrow(cluster_var_df))
  }else {
    stop("var_weight_type must be one of num_units, ess_units,
         uniform")
  }
    var_cluster <- cluster_var_df$var_cluster
    weighted_var = weighted.mean(var_cluster, w = weight)
  return( weighted_var )
}

#' Calculate the standard error estimate
#'
#' @param N_T number of treated
#' @param ESS_C effective size of controls
#' @param sigma_hat The estimated error standard deviation
#' @return The plug-in standard error estimate of the matching estimator
get_plug_in_SE <- function(N_T, ESS_C, sigma_hat) {
  sqrt((1 / N_T + 1 / ESS_C)) * sigma_hat
}

#' Main function: Estimate the variance from the plug-in estimator
#'
#' @param matches_table The data frame of the matched table
#' @param outcome Name of the outcome variable (default "Y")
#' @param treatment Name of the treatment variable (default "Z")
#' @return A tibble with SE, sigma_hat, N_T, and N_C_tilde
#' @export
get_se_AE_table <- function(
    matches_table,
    outcome = "Y",
    treatment = "Z",
    var_weight_type = "ess_units") {

  # Step 1: Filter the matched data to retain only control units with at least 2 in subclass
  matches_filtered <-
    matches_table %>%
    filter(!!sym(treatment) == 0) %>%
    group_by(subclass) %>%
    filter(n() >= 2) %>%
    ungroup()

  # Step 2: Calculate the variance for each subclass
  cluster_var_df <-
    calculate_subclass_variances(
      matches_filtered = matches_filtered,
      outcome = outcome
    )

  # Step 3: Calculate the weighted average variance
  weighted_var <-
    calculate_weighted_variance(
      cluster_var_df = cluster_var_df,
      var_weight_type = var_weight_type
    )
  sigma_hat <- sqrt(weighted_var)

  # Step 4: Calculate N_T and effective size of controls (N_C_tilde)
  Ns <- calc_N_T_N_C(matches_table)

  # Step 5: Calculate the plug-in standard error
  SE <- get_plug_in_SE(
    N_T = Ns$N_T,
    ESS_C = Ns$N_C_tilde,
    sigma_hat = sigma_hat
    )

  return(tibble(
    SE = SE,
    sigma_hat = sigma_hat,
    N_T = Ns$N_T,
    N_C_tilde = Ns$N_C_tilde
  ))
}

get_se_AE <- function(matches, outcome = "Y", treatment = "Z"){
  if ( is.csm_matches( matches ) ) {
    matches <- full_unit_table(matches)
  }

  get_se_AE_table(
    matches_table = matches,
    outcome = "Y",
    treatment = "Z"
  )
}

#' Estimate ATT
#'
#' Calculate ATT and associated standard error using the weighting
#' method from paper with locally estimated residual variation.
#'
#' @export
get_ATT_estimate <- function( scmatch, outcome = "Y" ) {

  ATE = get_att_ests( scmatch, outcome = outcome )
  se = get_se_AE( scmatch, outcome = outcome )
  se$ATE = ATE

  se %>% relocate( ATE ) %>%
    mutate(  t = ATE/SE )
}
