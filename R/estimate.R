# R/estimate.R



#' Estimate the ATT from a matched dataframe
#'
#' Given a matched dataset, calculate the estimated ATT
#'
#' @param csm A matched dataset, either the dataframe of treatment and
#'   control units, or a csm_matches object.  If dataframe, needs the
#'   treatment, outcome, and a "weights" column.
#' @param Either name of weight column or a numeric vector of weights
#'   can be provided.
#'
#' @noRd
get_att_point_est <- function(csm,
                              treatment = "Z",
                              outcome = "Y",
                              weights = "weights") {

  if ( is.csm_matches(csm) ) {
    csm <- result_table( csm, return="sc_units", outcome=outcome )
  }

  if ( is.numeric( weights ) ) {
    csm$weights <- weights
  } else if ( weights != "weights" ) {
    stopifnot( weights %in% names(csm) )
    csm$weights <- csm[[weights]]
  }

  stopifnot( all( c(treatment, outcome) %in% names(csm) ) )
  stopifnot( "weights" %in% names(csm) )

  csm$Z = csm[[treatment]]
  csm$Y = csm[[outcome]]


  csm %>%
    dplyr::group_by(Z) %>%
    dplyr::summarize(mn = sum(Y*weights) / sum(weights), .groups="drop") %>%
    dplyr::summarize(est = dplyr::last(mn) - dplyr::first(mn), .groups="drop") %>%
    dplyr::pull(est)
}


# get_att_ests_from_wt<-
#   function(df, input_wt)



calc_N_T_N_C <- function(preds_csm, treatment = "Z"){
  if ( is.csm_matches( preds_csm ) ) {
    treatment <- params(preds_csm)$treatment
    preds_csm <- result_table(preds_csm)
  }

  N_T <- nrow(preds_csm %>%
                dplyr::filter(.data[[treatment]] != 0))

  tmp <- preds_csm %>%
    dplyr::filter(.data[[treatment]]==0) %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(w_i = sum(weights), .groups="drop")

  N_C_tilde <- ess( tmp$w_i )


  return(list(N_T = N_T,
              N_C_tilde = N_C_tilde))
}


#' Effective Sample Size
#'
#' Calculate the effective sample size given a vector of unit weights
#'
#' @param weights A numeric vector of weights
#' @return The effective sample size
#' @export
#'
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


# #' Main function: Estimate the variance from the plug-in estimator
# #'
# #' @param matches_table The data frame of the matched table
# #' @param outcome Name of the outcome variable (default "Y")
# #' @param treatment Name of the treatment variable (default "Z")
# #' @param var_weight_type The way that cluster variances are averaged
# #' "num_units": weight by number of units in the subclass
# #' "ess_units": weight effective size of units in the subclass
# #' "uniform: weight each cluster equally
# #' @return A tibble with SE, sigma_hat, N_T, and N_C_tilde
# #' @export
# get_se_AE_table <- function(
    #     matches_table,
#   outcome = "Y",
#   treatment = "Z",
#   var_weight_type = "ess_units") {
#
# weighted_var <- get_pooled_variance(
#   matches_table = matches_table,
#   outcome = outcome,
#   treatment = treatment,
#   var_weight_type = var_weight_type)
#
# sigma_hat <- sqrt(weighted_var)
#
# # Step 4: Calculate N_T and effective size of controls (N_C_tilde)
# Ns <- calc_N_T_N_C(matches_table)
#
# # Step 5: Calculate the plug-in standard error
# SE <- get_plug_in_SE(
#   N_T = Ns$N_T,
#   ESS_C = Ns$N_C_tilde,
#   sigma_hat = sigma_hat
#   )
#




# get_se_AE <- function(matches,
#                       outcome = "Y",
#                       treatment = "Z",
#                       var_weight_type = "ess_units"){
#
#   if ( is.csm_matches( matches ) ) {
#     matches <- result_table(matches)
#   }
#
#   get_se_AE_table(
#     matches_table = matches,
#     outcome = outcome,
#     treatment = treatment,
#     var_weight_type = var_weight_type
#   )
# }



#' Get the standard error using the OR bootstrap approach
#'
#' This function estimates the standard error of the ATT using a bootstrap approach
#' based on residuals from the OR method.
#'
#' @param matches A CSM match object (R S3 object) or data frame
#' @param outcome Name of the outcome variable (default: "Y")
#' @param treatment Name of the treatment variable (default: "Z")
#' @param boot_mtd The bootstrap method to use. Options are "Bayesian", "wild", or "naive".
#' @param B Number of bootstrap samples (default: 250)
#' @param use_moving_block Whether to use a moving block bootstrap (default: FALSE)
#' @param seed_addition Additional seed value to ensure reproducibility (default: 11)
#' @param block_size Block size for bootstrap (default: NA, automatically chosen)
#'
#' @return A tibble with measurement error variance (V_E), standard error (SE), bootstrap confidence intervals, and sample sizes.
#' @export
#'
get_measurement_error_variance_OR <- function(matches,
                                              outcome = "Y",
                                              treatment = "Z",
                                              boot_mtd = "wild",
                                              B = 250,
                                              use_moving_block = FALSE,
                                              seed_addition = 11,
                                              block_size = NA) {
  # Convert to data frame if needed
  if (is.csm_matches(matches)) {
    matches <- result_table(matches)
  }

  # Check if bias column exists, if not set bias to 0
  if (!"bias" %in% colnames(matches)) {
    matches$bias <- 0
  }

  tmp0 <- matches %>%
    mutate(Y_bias_corrected = !!sym(outcome) - bias) %>%
    group_by(subclass, !!sym(treatment)) %>%
    summarize(mn = sum(Y_bias_corrected * weights), .groups = "drop")

  tmp <- tmp0 %>%
    group_by(subclass) %>%
    summarise(tilde_tau = last(mn) - first(mn), .groups = "drop")

  tilde_tau <- tmp$tilde_tau
  mean_tilde_tau <- mean(tilde_tau)
  tilde_tau_resids <- tilde_tau - mean_tilde_tau

  if (boot_mtd %in% c("Bayesian", "wild", "naive")) {
    boot_ci <- make_bootstrap_ci(
      boot_mtd,
      use_moving_block = use_moving_block
    )
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
  } else {
    stop("boot_mtd must be one of: 'Bayesian', 'wild', 'naive'")
  }

  # Calculate sample sizes
  Ns <- calc_N_T_N_C(matches, treatment = treatment)

  # Calculate V_E equivalent (variance per unit)
  V_E <- sd_boot^2

  return(tibble(
    V_E = V_E,
    SE = sd_boot,
    N_T = Ns$N_T,
    ESS_C = Ns$N_C_tilde,
    CI_lower = CI_lower,
    CI_upper = CI_upper
  ))
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
  sample_sizes <- calc_N_T_N_C(matches_table, treatment = treatment)
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


#' Estimate the measurement error variance component (V_E) under heterogeneous errors
#'
#' Calculates the measurement error variance component allowing for heterogeneous
#' variances across matched subclasses. This uses subclass-specific outcome variances
#' and propagates them through the weighting structure to estimate V_E.
#'
#' @param matches_table The data frame of the matched table (e.g., output from `full_unit_table`)
#' @param outcome Name of the outcome variable (default "Y")
#' @param treatment Name of the treatment variable (default "Z")
#'
#' @return A tibble with measurement error variance (V_E), pooled sigma_hat, N_T, and ESS_C
#' @export
get_measurement_error_variance_het <- function(
    matches_table,
    outcome = "Y",
    treatment = "Z",
    cluster_comb_mtd = "sample") {

  # Keep only subclasses with at least 3 units to estimate within-subclass variance
  matches_filtered <-
    matches_table %>%
    group_by(subclass) %>%
    filter(n() >= 3) %>%
    ungroup()

  # Subset to control units (Z == 0) for estimating subclass variances
  matches_filtered_co <-
    matches_filtered %>%
    filter(!!sym(treatment) == 0)

  # Estimate outcome variance within each control subclass
  cluster_var_df <-
    calculate_subclass_variances(
      matches_filtered = matches_filtered_co,
      outcome = outcome
    )

  # Attach subclass variances to all units in the matched table
  matches_with_var <-
    matches_filtered %>%
    left_join(cluster_var_df, by = "subclass")

  # Compute total weight and average subclass variance per unit
  if (nrow(matches_with_var) == 0) {
    # Create an empty tibble with the correct columns
    var_calc_df <- tibble(
      id = integer(0),
      !!sym(treatment) := integer(0),
      total_wt = numeric(0),
      total_wt_squared = numeric(0),
      avg_var_cluster = numeric(0),
      rand_var_cluster = numeric(0)
    )
  } else {
    # Compute total weight and average subclass variance per unit
    var_calc_df <- matches_with_var %>%
      group_by(id, !!sym(treatment)) %>%
      summarise(
        total_wt = sum(weights),
        total_wt_squared = sum(weights^2),
        # Add na.rm = TRUE for cases like 'mock_matches_one_control'
        avg_var_cluster = mean(var_cluster, na.rm = TRUE),
        # Add a safety check for sample()
        rand_var_cluster = if(length(var_cluster) > 0 && !all(is.na(var_cluster))) {
          sample(var_cluster[!is.na(var_cluster)], 1)
        } else {
          NA_real_
        },
        .groups = 'drop'
      )
  }

  # Calculate sample sizes
  sample_sizes <- calc_N_T_N_C(matches_table)
  N_T <- sample_sizes$N_T
  ESS_C <- sample_sizes$N_C_tilde

  # Compute V_E accounting for heteroskedasticity
  if (cluster_comb_mtd == "sample") {
    V_E_het <- sum(var_calc_df$total_wt^2 * var_calc_df$rand_var_cluster) / N_T^2
  } else if (cluster_comb_mtd == "average") {
    V_E_het <- sum(var_calc_df$total_wt^2 * var_calc_df$avg_var_cluster) / N_T^2
  } else {
    stop("Invalid value for `cluster_comb_mtd`. Use either 'sample' or 'average'.")
  }

  # Also return pooled variance estimate for completion
  weighted_var <- get_pooled_variance(
    matches_table = matches_table,
    outcome = outcome,
    treatment = treatment)

  sigma_hat <- sqrt(weighted_var)

  return(list(
    V_E = V_E_het,
    sigma_hat = sigma_hat,
    N_T = N_T,
    ESS_C = ESS_C,
    var_calc_df = var_calc_df
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
#' @param variance_method Method for calculating measurement error variance:
#'   "pooled": use get_measurement_error_variance (default)
#'   "bootstrap": use get_measurement_error_variance_OR
#' @param boot_mtd Bootstrap method when variance_method = "bootstrap" (default: "wild")
#' @param B Number of bootstrap samples when variance_method = "bootstrap" (default: 250)
#' @param seed_addition Additional seed for bootstrap (default: 11)
#' @return A tibble with total variance (V), measurement error variance (V_E),
#'   population heterogeneity variance (V_P), and other relevant statistics
#' @export
get_total_variance <- function(
    matches,
    outcome = "Y",
    treatment = "Z",
    var_weight_type = "ess_units",
    variance_method = "pooled",
    boot_mtd = "wild",
    B = 250,
    seed_addition = 11,
    cluster_comb_mtd = "average",
    df = NULL,
    ...) {

  # Convert to data frame if needed
  if (is.csm_matches(matches)) {
    matches_df <- full_unit_table(matches)
  } else {
    matches_df <- matches
  }

  # --- Capture extra arguments ---
  extra_args <- list(...)

  # Get measurement error variance estimate (V_E) using specified method
  if (variance_method == "pooled") {
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
  } else if (variance_method == "bootstrap") {
    v_e_result <- get_measurement_error_variance_OR(
      matches = matches_df,
      outcome = outcome,
      treatment = treatment,
      boot_mtd = boot_mtd,
      B = B,
      seed_addition = seed_addition
    )
    V_E <- v_e_result$V_E
    # For bootstrap method, we don't have sigma_hat, so we approximate
    # Essentially we are using V_E as the whole V for bootstrap, which is desired.
    sigma_hat_squared <- V_E * (v_e_result$N_T + v_e_result$ESS_C) /
      (1/v_e_result$N_T + 1/v_e_result$ESS_C)
    N_T <- v_e_result$N_T
    ESS_C <- v_e_result$ESS_C
  } else if (variance_method == "pooled_het"){
    v_e_result <- get_measurement_error_variance_het(
      matches_table = matches_df,
      outcome = outcome,
      treatment = treatment,
      cluster_comb_mtd = cluster_comb_mtd
    )
    V_E <- v_e_result$V_E
    sigma_hat_squared <- v_e_result$sigma_hat^2
    N_T <- v_e_result$N_T
    ESS_C <- v_e_result$ESS_C
  } else if (variance_method == "ai06") {
    if (is.null(df)) {
      stop("The 'ai06' variance_method requires the full 'df' argument.")
    }
    if (is.null(extra_args$M)) {
      stop("The 'ai06' variance_method requires the 'M' argument.")
    }

    v_e_result_list <- get_variance_AI06(
      df = df,
      scmatch = matches, # Pass the original S3 object
      outcome = outcome,
      treatment = treatment,
      ... # Pass M, covs, scaling, metric, etc.
    )
    # Unpack complex object
    v_e_result <- v_e_result_list %>% select(-var_calc_df)


    V_E <- v_e_result_list$V_E
    sigma_hat_squared <- v_e_result_list$sigma_hat^2
    N_T <- v_e_result_list$N_T
    ESS_C <- v_e_result_list$ESS_C
    var_calc_df <- v_e_result_list$var_calc_df[[1]]
    # var_calc_df <- v_e_result_list$var_calc_df # Extract the data frame
    v_e_result <- v_e_result_list # for sigma_hat lookup later

  } else {
    # Update error message
    stop("variance_method must be 'pooled', 'bootstrap', 'pooled_het', or 'ai06'")
  }

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

  squared_deviations_se_ver <- var(individual_effects$indiv_effect)


  # Calculate the correction term for control unit weights
  if (variance_method == "pooled" | variance_method == "bootstrap"){
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

    V_correction <- sigma_hat_squared * control_weight_correction
  } else if (variance_method == "pooled_het") {
    var_calc_df <- v_e_result$var_calc_df


    if (cluster_comb_mtd == "average"){
      V_correction <- var_calc_df %>%
        filter(!!sym(treatment) == 0) %>%
        summarize(
          V_correction_term = sum((total_wt^2 - total_wt_squared) * avg_var_cluster / N_T )
        ) %>%
        pull(V_correction_term)
    } else if (cluster_comb_mtd == "sample"){
      V_correction <- var_calc_df %>%
        filter(!!sym(treatment) == 0) %>%
        summarize(
          V_correction_term = sum((total_wt^2 - total_wt_squared) * sample_var_cluster / N_T )
        ) %>%
        pull(V_correction_term)
    } else {
      stop("Invalid value for `cluster_comb_mtd`. Use either 'sample' or 'average'.")
    }

  } else if (variance_method == "ai06") {
    # V_correction logic for ai06
    # Note: 'avg_var_cluster' is the column name we assigned to sigma_j^2
    if (is.null(var_calc_df)) stop("var_calc_df is NULL for ai06")

    V_correction <- as.tibble(var_calc_df) %>%
      filter(!!sym(treatment) == 0) %>%
      summarize(
        V_correction_term = sum((total_wt^2 - total_wt_squared) * avg_var_cluster / N_T, na.rm = TRUE)
      ) %>%
      pull(V_correction_term)

  } else {
    # Should be unreachable
    V_correction <- 0
  }


  # Calculate the total variance (V)
  V <- squared_deviations + V_correction

  # Calculate population heterogeneity variance (V_P)
  N_T_V_P <- V - N_T * V_E
  if (variance_method == "bootstrap"){
    V = N_T * V_E
    N_T_V_P = 0
  }



  SE <- sqrt(V) * 1 / sqrt(N_T)

  result <- tibble(
    V = V,
    V_E = V_E,
    V_P = N_T_V_P / N_T,
    SE = SE,
    N_T = N_T,
    ESS_C = ESS_C,
    squared_deviations = squared_deviations,
    squared_deviations_se_ver,
    V_correction = V_correction
  )

  # Add sigma_hat if available (from pooled method we have, from OR method we impute)
  if (variance_method == "pooled" | variance_method == "pooled_het" | variance_method == "ai06") {
    result$sigma_hat <- v_e_result$sigma_hat
  }else if (variance_method == "bootstrap") {
    result$sigma_hat <- sqrt(sigma_hat_squared)
  } else {
    result$sigma_hat <- NA
  }

  return(result)
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
#' @param variance_method Method for calculating measurement error variance:
#'   "pooled": use get_measurement_error_variance (default)
#'   "bootstrap": use get_measurement_error_variance_OR
#' @param boot_mtd Bootstrap method when variance_method = "bootstrap" (default: "wild")
#' @param B Number of bootstrap samples when variance_method = "bootstrap" (default: 250)
#' @param seed_addition Additional seed for bootstrap (default: 11)
#' @param df The *full* original data frame. Required for 'ai06' method.
#' @param ... Additional arguments passed to get_total_variance.
#' @return A tibble with ATT estimate, standard error, t-statistic, and variance components
#' @export
get_ATT_estimate <- function(
    scmatch,
    treatment = "Z",
    outcome = "Y",
    var_weight_type = "ess_units",
    variance_method = "pooled",
    boot_mtd = "wild",
    B = 250,
    seed_addition = 11,
    cluster_comb_mtd = "average",
    df = NULL, # Add df
    ...) {

  # Calculate ATT point estimate
  ATT <- get_att_point_est(scmatch,
                           treatment = treatment,
                           outcome = outcome)

  # Calculate total variance and standard error
  variance_results <- get_total_variance(
    matches = scmatch,
    treatment = treatment,
    outcome = outcome,
    var_weight_type = var_weight_type,
    variance_method = variance_method,
    boot_mtd = boot_mtd,
    B = B,
    seed_addition = seed_addition,
    cluster_comb_mtd = cluster_comb_mtd,
    df = df, # Pass df
    ... # Pass ... (e.g., M)
  )

  # Add ATT to results and calculate t-statistic
  variance_results %>%
    mutate(ATT = ATT) %>%
    relocate(ATT) %>%
    mutate(t = ATT/SE)
}


#' Estimate ATT and SE
#'
#' Calculate ATT and associated standard error using the total
#' variance estimator from the paper that accounts for both
#' measurement error and population heterogeneity.
#'
#' @param scmatch The CSM match object, an R S3 object, or a data
#'   frame of the full results.
#' @param treatment Name of the treatment variable (default "Z")
#' @param outcome Name of the outcome variable (default "Y")
#' @param var_weight_type The way that cluster variances are averaged
#'   (default "ess_units"): "num_units": weight by number of units in
#'   the subclass "ess_units": weight by effective sample size of
#'   units in the subclass "uniform": weight each cluster equally
#' @param variance_method Method for calculating measurement error
#'   variance: "pooled": use get_measurement_error_variance (default)
#'   "bootstrap": use get_measurement_error_variance_OR
#' @param boot_mtd Bootstrap method when variance_method = "bootstrap"
#'   (default: "wild")
#' @param B Number of bootstrap samples when variance_method =
#'   "bootstrap" (default: 250)
#' @param seed_addition Additional seed for bootstrap (default: 11)
#' @param feasible_only Logical indicating whether to use only feasible
#'   matches (default: FALSE)
#' @export
#'
#' @return A tibble with ATT estimate (ATT), standard error (SE),
#'   total variance (V), measurement error variance (V_E), population
#'   heterogeneity variance (V_P), and other relevant statistics
#'
#' @export
#'
estimate_ATT <- function(
    matches,
    outcome = NULL,
    treatment = "Z",
    var_weight_type = "ess_units",
    variance_method = "pooled",
    boot_mtd = "wild",
    B = 250,
    seed_addition = 11,
    feasible_only = FALSE ) {

  # Convert to data frame if needed
  if (is.csm_matches(matches)) {
    matches_df <- result_table(matches, feasible_only = feasible_only )
    treatment = params(matches)$treatment
  } else {
    matches_df <- matches
    if ( feasible_only ) {
      stopifnot( "feasible" %in% colnames( matches ) )
      matches_df <- matches %>%
        filter( feasible == 1 )
    }
  }

  # return blank results df if no units matched
  if ( nrow(matches_df) == 0 ) {
    result <- tibble(
      ATT = NA,
      SE = NA,
      N_T = NA,
      ESS_C = NA,
      t = NA,
      V = NA,
      V_E = NA,
      V_P = NA,
      sigma_hat = NA
    )
    return(result)
  }


  if ( n_distinct( matches_df[[treatment]] ) != 2 ) {
    stop(glue::glue( "treatment variable `{treatment}` must be binary (0/1)") )
  }

  if ( is.null(outcome) ) {
    if ( variance_method != "pooled" ) {
      stop( glue::glue( "outcome must be specified when using {variance_method} variance method") )
    }

    matches_df$Y = 0
    v_e_result <- get_measurement_error_variance(
      matches_table = matches_df,
      outcome = "Y",
      treatment = treatment,
      var_weight_type = var_weight_type
    )
    tt <- tibble( N_T = v_e_result$N_T,
                  ESS_C = v_e_result$ESS_C )
    return( tt )
  }


  # Get measurement error variance estimate (V_E) using specified method
  if (variance_method == "pooled") {
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
  } else if (variance_method == "bootstrap") {
    v_e_result <- get_measurement_error_variance_OR(
      matches = matches_df,
      outcome = outcome,
      treatment = treatment,
      boot_mtd = boot_mtd,
      B = B,
      seed_addition = seed_addition
    )
    V_E <- v_e_result$V_E
    # For bootstrap method, we don't have sigma_hat, so we approximate
    sigma_hat_squared <- V_E * (v_e_result$N_T + v_e_result$ESS_C) /
      (1/v_e_result$N_T + 1/v_e_result$ESS_C)
    N_T <- v_e_result$N_T
    ESS_C <- v_e_result$ESS_C
  } else {
    stop("variance_method must be either 'pooled' or 'bootstrap'")
  }

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

  # Calculate the total variance (V)
  V <- squared_deviations + sigma_hat_squared * control_weight_correction

  # Calculate population heterogeneity variance (V_P)
  V_P <- V - N_T * V_E
  if (variance_method == "bootstrap"){
    V = N_T * V_E
    V_P = 0
  }

  SE <- sqrt(V) * 1 / sqrt(N_T)

  result <- tibble(
    ATT = att_estimate,
    SE = SE,
    N_T = N_T,
    ESS_C = ESS_C,
    t = att_estimate / SE,
    V = V,
    V_E = V_E,
    V_P = V_P
  )

  # Add sigma_hat if available (from pooled method we have, from OR method we impute)
  if (variance_method == "pooled") {
    result$sigma_hat <- v_e_result$sigma_hat
  }else if (variance_method == "bootstrap") {
    result$sigma_hat <- sqrt(sigma_hat_squared)
  } else {
    result$sigma_hat <- NA
  }

  return(result)
}


