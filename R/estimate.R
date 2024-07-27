


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

get_se_AE <- function(preds_csm, outcome = "Y"){
  # 1. Get debiased units; Get the subclasses

  # TODO: Fix this.  where is hat_mu_0 set?
  #preds_csm <- preds_csm$result %>%
  #  mutate(Y_bias_corrected = Y - hat_mu_0)

  # 2. Filter the controls
  #     and the subclasses with n_controls >= 2
  preds_csm_filtered <-
    full_unit_table(preds_csm) %>%
    filter(Z == 0) %>%
    group_by(subclass) %>%
    filter(n() >= 2) %>%
    ungroup()

  # 3. For each filtered subclass,
  #   calculate the cluster residual s_j^2 = se(debiased_units)
  weighted_var <- function(x, wt) {
    n <- length(x)
    wm <- weighted.mean(x, wt)
    sum(wt * (x - wm)^2) * n / (n-1)
  }

  weighted_se <- function(x, wt) {
    sqrt(weighted_var(x, wt) / length(x))
  }

  cluster_var_df <-
    preds_csm_filtered %>%
    group_by(subclass) %>%
    summarise(nj = n(),
              w_nj = ess(weights),
              var_cluster = var(Y), .groups="drop")
  # var_cluster = weighted_var(Y_bias_corrected, weights))

  # 4. Get the weighted average of s_j^2, weighted by n_j, number of
  # units in the subclass
  #
  # NOTE: Could weight by nj or w_nj.  w_nj takes into account how
  # much uncertainty is in each group, and thus might perform better
  # with heteroskedasticity?
  weighted_var_df <- cluster_var_df %>%
    summarise(weighted_var = weighted.mean(var_cluster, w = w_nj), .groups="drop")
  sigma_hat <- sqrt(weighted_var_df$weighted_var)

  # 5. calculate N_T and N_C
  Ns <- calc_N_T_N_C(preds_csm)

  # 6. Calculate the variance of the estimator
  res <- sqrt((1/Ns$N_T + 1/Ns$N_C_tilde)) * sigma_hat
  return( tibble( SE = res,
                  sigma_hat = sigma_hat,
                  N_T = Ns$N_T,
                  N_C_tilde = Ns$N_C_tilde) )
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
