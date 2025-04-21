


# functions for estimating effects
# The following two functions should be combined
#   once we get to know the input of get_att_point_est
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
#' @param matched_df A matched dataset, either the dataframe of
#'   treatment and control units, or a csm_matches object.  If
#'   dataframe, needs the treatment, outcome, and a "weights" column.
#'
#' @export
get_att_point_est <- function(matched_df, treatment = "Z", outcome = "Y") {

  if ( is.csm_matches(matched_df) ) {
    matched_df <- result_table( matched_df, "sc_units" )
  }
  stopifnot( all( c(treatment, outcome) %in% names(matched_df) ) )
  stopifnot( "weights" %in% names(matched_df) )

  matched_df$Z = matched_df[[treatment]]
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
#' @param var_weight_type The way that cluster variances are averaged
#' "num_units": weight by number of units in the subclass
#' "ess_units": weight effective size of units in the subclass
#' "uniform: weight each cluster equally
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


get_pooled_variance <- function(
    matches_table,
    outcome = "Y",
    treatment = "Z",
    var_weight_type = "ess_units"){
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
  return(weighted_var)

}


#' #' Main function: Estimate the variance from the plug-in estimator
#' #'
#' #' @param matches_table The data frame of the matched table
#' #' @param outcome Name of the outcome variable (default "Y")
#' #' @param treatment Name of the treatment variable (default "Z")
#' #' @param var_weight_type The way that cluster variances are averaged
#' #' "num_units": weight by number of units in the subclass
#' #' "ess_units": weight effective size of units in the subclass
#' #' "uniform: weight each cluster equally
#' #' @return A tibble with SE, sigma_hat, N_T, and N_C_tilde
#' #' @export
#' get_se_AE_table <- function(
#'     matches_table,
#'     outcome = "Y",
#'     treatment = "Z",
#'     var_weight_type = "ess_units") {
#'
#'   weighted_var <- get_pooled_variance(
#'     matches_table = matches_table,
#'     outcome = outcome,
#'     treatment = treatment,
#'     var_weight_type = var_weight_type)
#'
#'   sigma_hat <- sqrt(weighted_var)
#'
#'   # Step 4: Calculate N_T and effective size of controls (N_C_tilde)
#'   Ns <- calc_N_T_N_C(matches_table)
#'
#'   # Step 5: Calculate the plug-in standard error
#'   SE <- get_plug_in_SE(
#'     N_T = Ns$N_T,
#'     ESS_C = Ns$N_C_tilde,
#'     sigma_hat = sigma_hat
#'     )
#'
#'   return(tibble(
#'     SE = SE,
#'     sigma_hat = sigma_hat,
#'     N_T = Ns$N_T,
#'     N_C_tilde = Ns$N_C_tilde
#'   ))
#' }
#'
#' #' Get the standard error using the weighting approach
#' #'
#' #' Method taken from the balancing weights literature.
#' #'
#' #' @param matches The CSM match object, an R S3 object
#' #' @param outcome Name of the outcome variable (default "Y")
#' #' @param treatment Name of the treatment variable (default "Z")
#' #' @param var_weight_type The way that cluster variances are averaged
#' #'   "num_units": weight by number of units in the subclass
#' #'   "ess_units": weight effective size of units in the subclass
#' #'   "uniform: weight each cluster equally
#' #'
#' #' @return A tibble with SE, sigma_hat, N_T, and N_C_tilde
#' #' @export
#' #'
#' get_se_AE <- function(matches,
#'                       outcome = "Y",
#'                       treatment = "Z",
#'                       var_weight_type = "ess_units"){
#'   if ( is.csm_matches( matches ) ) {
#'     matches <- full_unit_table(matches)
#'   }
#'
#'   get_se_AE_table(
#'     matches_table = matches,
#'     outcome = outcome,
#'     treatment = treatment,
#'     var_weight_type = var_weight_type
#'   )
#' }

get_se_OR <- function(matches,
                      outcome = "Y",
                      treatment = "Z",
                      boot_mtd = "wild",
                      B = 250,
                      use_moving_block=F,
                      seed_addition = 11,
                      block_size = NA){
  # matches <- mock_matches; outcome = "Y"; treatment = "Z"; boot_mtd = "wild"
  if ( is.csm_matches( matches ) ) {
    matches <- full_unit_table(matches)
  }
  tmp0 <- matches %>%
    mutate(Y_bias_corrected = Y) %>%
    group_by(subclass, Z) %>%
    summarize(mn = sum(Y_bias_corrected * weights), .groups = "drop")

  tmp <- tmp0 %>%
    group_by(subclass) %>%
    summarise(tilde_tau = last(mn) - first(mn), .groups = "drop")

  tilde_tau <- tmp$tilde_tau
  mean_tilde_tau <- mean(tilde_tau)
  tilde_tau_resids <- tilde_tau - mean_tilde_tau

  if (boot_mtd == "Bayesian" || boot_mtd == "wild" || boot_mtd == "naive-resid"){
    boot_ci <- make_bootstrap_ci(
      boot_mtd,
      use_moving_block=use_moving_block
    )  # or "wild" or "naive-resid"
    results <- boot_ci(
      resids = tilde_tau_resids,
      mean_est = mean_tilde_tau,
      B = B,
      seed_addition = seed_addition,
      block_size = block_size
    )

    CI_lower <- results$ci_lower
    CI_upper <- results$ci_upper
    sd_boot <- results$sd
  }
  Ns <- calc_N_T_N_C(matches)

  return( tibble(
    SE = sd_boot,
    sigma_hat = NA,
    N_T = Ns$N_T,
    N_C_tilde = Ns$N_C_tilde,
    CI_lower = CI_lower,
    CI_upper = CI_upper
  ) )
}


#' Estimate the measurement error variance component (V_E)
#'
#' Calculates the measurement error variance component using the pooled variance
#' approach described in the paper. This estimates the variance component
#' due to the noise in outcomes.
#'
#' @param matches_table The data frame of the matched table
#' @param outcome Name of the outcome variable (default "Y")
#' @param treatment Name of the treatment variable (default "Z")
#' @param var_weight_type The way that cluster variances are averaged:
#'   "num_units": weight by number of units in the subclass
#'   "ess_units": weight by effective sample size of units in the subclass
#'   "uniform": weight each cluster equally
#' @return A tibble with measurement error variance (V_E), sigma_hat, N_T, and ESS_C
#' @export
get_measurement_error_variance <- function(
    matches_table,
    outcome = "Y",
    treatment = "Z",
    var_weight_type = "ess_units") {

  # Get pooled variance estimate (S^2 in Equation 10)
  weighted_var <- get_pooled_variance(
    matches_table = matches_table,
    outcome = outcome,
    treatment = treatment,
    var_weight_type = var_weight_type)

  sigma_hat <- sqrt(weighted_var)

  # Calculate N_T and effective sample size of controls (ESS_C)
  sample_sizes <- calc_N_T_N_C(matches_table)
  N_T <- sample_sizes$N_T
  ESS_C <- sample_sizes$N_C_tilde

  # Calculate the measurement error variance V_E (Equation 9)
  V_E <- weighted_var * (1/N_T + 1/ESS_C)

  return(tibble(
    V_E = V_E,
    sigma_hat = sigma_hat,
    N_T = N_T,
    ESS_C = ESS_C
  ))
}

#' Calculate the total variance estimator (V)
#'
#' Implements the total variance estimator from the paper, which accounts for
#' both measurement error variance (V_E) and population heterogeneity variance (V_P).
#'
#' @param matches The CSM match object, an R S3 object
#' @param outcome Name of the outcome variable (default "Y")
#' @param treatment Name of the treatment variable (default "Z")
#' @param var_weight_type The way that cluster variances are averaged:
#'   "num_units": weight by number of units in the subclass
#'   "ess_units": weight by effective sample size of units in the subclass
#'   "uniform": weight each cluster equally
#' @return A tibble with total variance (V), measurement error variance (V_E),
#'   population heterogeneity variance (V_P), and other relevant statistics
#' @export
get_total_variance <- function(
    matches,
    outcome = "Y",
    treatment = "Z",
    var_weight_type = "ess_units") {

  # Convert to data frame if needed
  if (is.csm_matches(matches)) {
    matches_df <- full_unit_table(matches)
  } else {
    matches_df <- matches
  }

  # Get measurement error variance estimate (V_E)
  v_e_result <- get_measurement_error_variance(
    matches_table = matches_df,
    outcome = outcome,
    treatment = treatment,
    var_weight_type = var_weight_type
  )

  V_E <- v_e_result$V_E
  sigma_hat_squared <- v_e_result$sigma_hat^2
  N_T <- v_e_result$N_T
  ESS_C <- v_e_result$ESS_C

  # Calculate treatment effect estimate (tau_hat)
  att_estimate <- get_att_point_est(
    matches_df,
    treatment = treatment,
    outcome = outcome
  )

  # Calculate individual treatment effects
  # First, get the predicted counterfactual outcomes for each treated unit
  tx_units <- matches_df %>%
    filter(!!sym(treatment) == 1) %>%
    select(id, subclass, !!sym(outcome)) %>%
    rename(Y_t = !!sym(outcome))

  # Get the weighted control outcomes for each subclass (Y_hat(0))
  control_outcomes <- matches_df %>%
    filter(!!sym(treatment) == 0) %>%
    group_by(subclass) %>%
    summarize(
      Y_hat_0 = sum(!!sym(outcome) * weights) / sum(weights),
      .groups = "drop"
    )

  # Join to get individual treatment effects
  individual_effects <- tx_units %>%
    left_join(control_outcomes, by = "subclass") %>%
    mutate(indiv_effect = Y_t - Y_hat_0)

  # Calculate the first component of the total variance estimator
  # (Empirical squared deviations)
  squared_deviations <- sum((individual_effects$indiv_effect - att_estimate)^2) / N_T
  # print("Individual treatment effects:" )
  # print(individual_effects$indiv_effect)
  # print("Average treatment effect:" )
  # print(mean(individual_effects$indiv_effect))
  # print("Squared deviations:" )
  # print(squared_deviations)

  # Calculate the correction term for control unit weights
  control_weight_correction <- matches_df %>%
    filter(!!sym(treatment) == 0) %>%
    group_by(id) %>%
    summarize(
      sum_weights = sum(weights),
      sum_squared_weights = sum(weights^2),
      .groups = "drop"
    ) %>%
    summarize(
      correction_term = sum((sum_weights^2 - sum_squared_weights)) / N_T
    ) %>%
    pull(correction_term)
  # print("Control weight correction term:" )
  # print(control_weight_correction)
  # Calculate the total variance (V)
  V <- squared_deviations + sigma_hat_squared * control_weight_correction

  # Calculate population heterogeneity variance (V_P)
  V_P <- V - N_T * V_E

  SE <- sqrt(V) * 1 / sqrt(N_T)

  return(tibble(
    V = V,
    V_E = V_E,
    V_P = V_P,
    SE = SE,
    sigma_hat = v_e_result$sigma_hat,
    N_T = N_T,
    ESS_C = ESS_C
  ))
}

#' Estimate ATT with correct variance estimation
#'
#' Calculate ATT and associated standard error using the total variance estimator
#' from the paper that accounts for both measurement error and population heterogeneity.
#'
#' @param scmatch The CSM match object, an R S3 object
#' @param treatment Name of the treatment variable (default "Z")
#' @param outcome Name of the outcome variable (default "Y")
#' @param var_weight_type The way that cluster variances are averaged (default "ess_units")
#' @return A tibble with ATT estimate, standard error, t-statistic, and variance components
#' @export
get_ATT_estimate <- function(
    scmatch,
    treatment = "Z",
    outcome = "Y",
    var_weight_type = "ess_units") {

  # Calculate ATT point estimate
  ATT <- get_att_point_est(scmatch,
                           treatment = treatment,
                           outcome = outcome)

  # Calculate total variance and standard error
  variance_results <- get_total_variance(
    matches = scmatch,
    treatment = treatment,
    outcome = outcome,
    var_weight_type = var_weight_type
  )

  # Add ATT to results and calculate t-statistic
  variance_results %>%
    mutate(ATT = ATT) %>%
    relocate(ATT) %>%
    mutate(t = ATT/SE)
}
