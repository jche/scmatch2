

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
              N_C = sum(tmp$w_i > 0.00000001),
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
#' This should be removed, I think.
#'
#' @param N_T number of treated
#' @param ESS_C effective size of controls
#' @param sigma_hat The estimated error standard deviation
#' @return The plug-in standard error estimate of the matching estimator
#'
#' @noRd
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
#'   "pooled": use get_measurement_error_variance (default).
#'   "bootstrap": use get_measurement_error_variance_OR.
#'   "pooled_het": allow for heterogeneity.
#'   "ai06": the estimator of Abadie and Imbens.
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
    variance_method = c( "pooled", "pooled_het", "bootstrap", "ai06" ),
    boot_mtd = "wild",
    B = 250,
    seed_addition = 11,
    cluster_comb_mtd = "average",
    df = NULL,
    ...) {

  variance_method = match.arg(variance_method)

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
    v_e_result <- v_e_result_list %>%
      dplyr::select(-var_calc_df)


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
    dplyr::select(id, subclass, !!sym(outcome)) %>%
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




#' Estimate ATT and SE
#'
#' Calculate ATT and associated standard error using the specified
#' variance estimator.
#'
#' This will implement both finite-sample and superpopulation variance
#' estimation.  For superpopulation inference, it can also do the OR
#' bootstrap, if specified.
#'
#' @param matches The CSM match object (csm_matches), or a data frame
#'   of matched units (must contain treatment, outcome, and weights
#'   columns).
#' @param outcome Name of the outcome variable (default "Y").
#' @param treatment Name of the treatment variable (default "Z").
#' @param superpopulation If TRUE (default), use the superpopulation
#'   variance estimator via \code{get_total_variance()}.  If FALSE,
#'   use the finite-sample estimator via \code{get_finite_variance()}.
#' @param var_weight_type How cluster variances are averaged (default
#'   "ess_units"): "num_units" weights by number of units in the
#'   subclass; "ess_units" weights by effective sample size; "uniform"
#'   weights each cluster equally.
#' @param variance_method Variance method passed to
#'   \code{get_total_variance()}: "pooled" (default), "pooled_het",
#'   "bootstrap", or "ai06".  Ignored when \code{superpopulation =
#'   FALSE}.
#' @param boot_mtd Bootstrap method when variance_method = "bootstrap"
#'   (default: "wild").
#' @param B Number of bootstrap samples (default: 250).
#' @param seed_addition Additional seed for bootstrap (default: 11).
#' @param cluster_comb_mtd How per-unit variances are combined when
#'   variance_method = "pooled_het" (default: "average").
#' @param feasible_only Logical; use only feasible matches (default
#'   FALSE).
#' @param df The *full* original data frame. Required for 'ai06'
#'   method and for \code{get_finite_variance()} when
#'   \code{use_common_variance = FALSE}.
#' @param ... Additional arguments passed to the chosen variance
#'   function. For \code{get_total_variance()}: e.g. \code{M} for
#'   'ai06'. For \code{get_finite_variance()}: e.g.
#'   \code{use_common_variance}, \code{K}, \code{covs},
#'   \code{scaling}.
#' @return A tibble with ATT estimate (ATT), standard error (SE),
#'   t-statistic (t), total variance (V), measurement error variance
#'   (V_E), population heterogeneity variance (V_P), and other
#'   relevant statistics.
#'
#' @export
#'
estimate_ATT <- function(
    matches,
    outcome = "Y",
    treatment = "Z",
    superpopulation = FALSE,
    homoskedastic = FALSE,
    use_common_variance = TRUE,
    var_weight_type = "ess_units",
    variance_method = "pooled",
    boot_mtd = "wild",
    B = 250,
    seed_addition = 11,
    cluster_comb_mtd = "average",
    feasible_only = FALSE,
    df = NULL,
    ... ) {

  # Unpack CSM object or plain matched data frame
  if (is.csm_matches(matches)) {
    matches_df <- result_table(matches, feasible_only = feasible_only)
    treatment  <- params(matches)$treatment
  } else {
    matches_df <- matches
    if (feasible_only) {
      stopifnot("feasible" %in% colnames(matches))
      matches_df <- matches %>% filter(feasible == 1)
    }
  }

  # Return blank result if nothing matched
  if (nrow(matches_df) == 0) {
    return(tibble(
      ATT = NA_real_, SE = NA_real_, t = NA_real_,
      V = NA_real_, V_E = NA_real_, V_P = NA_real_,
      N_T = NA_real_, ESS_C = NA_real_, N_C = NA_real_, sigma_hat = NA_real_
    ))
  }

  if ( !( outcome %in% colnames(matches_df) ) ) {
    stop( glue::glue( "Outcome variable {outcome} not found in data" ) )
  }


  # Point estimate
  ATT <- get_att_point_est(matches_df,
                           treatment = treatment,
                           outcome   = outcome)


  # Variance + SE — two paths
  if (superpopulation) {
    variance_results <- get_total_variance(
      matches          = matches_df,
      outcome          = outcome,
      treatment        = treatment,
      var_weight_type  = var_weight_type,
      variance_method  = variance_method,
      boot_mtd         = boot_mtd,
      B                = B,
      seed_addition    = seed_addition,
      cluster_comb_mtd = cluster_comb_mtd,
      df               = df,
      ...
    )
  } else {
    if ( is.csm_matches( matches ) ) {
      df = full_unit_table(matches)
      pp = params(matches)
      scaling = pp$scaling
      covs = pp$covs
      metric = pp$metric
      id_name = pp$id_name
    }
    variance_results <- get_finite_variance(
      matches   = matches_df,
      outcome   = outcome,
      treatment = treatment,
      df        = df,
      homoskedastic = homoskedastic,
      use_common_variance = use_common_variance,
      scaling = scaling,
      covs = covs,
      metric = metric,
      id_name = id_name,
      ...
    )
  }

  # Attach ATT, t-statistic, and unique control count
  variance_results %>%
    mutate(ATT = ATT) %>%
    relocate(ATT) %>%
    mutate(t = ATT / SE)
}






#' Compute per-subclass s_t^2 and per-control-unit s_j^2
#'
#' For each treated unit's matched set (subclass), computes the variance of
#' control outcomes s_t^2.  Then, for each control unit j, averages the s_t^2
#' values over all subclasses that j appears in to obtain s_j^2.
#'
#' @param matches_table Full matched table (output of full_unit_table or result_table)
#' @param outcome Name of outcome variable (default "Y")
#' @param treatment Name of treatment variable (default "Z")
#' @return A list with:
#'   \describe{
#'     \item{s_t_sq}{Treatment unit variances. Data frame with columns \code{subclass} and \code{s_t_sq}.
#'       Only subclasses with at least 2 controls are included.}
#'     \item{s_j_sq}{Control unit variances. Data frame with columns \code{id} and \code{s_j_sq}.
#'       One row per control unit.  Units whose subclasses all have fewer than 2
#'       controls will have \code{NaN}.}
#'   }
#' @noRd
#'
calculate_s_j_sq <- function(matches_table, outcome = "Y", treatment = "Z") {
  matches_table <- matches_table %>%
    mutate(id = as.character(id))

  s_t_sq_df <- matches_table %>%
    filter(!!sym(treatment) == 0) %>%
    group_by(subclass) %>%
    filter(n() >= 2) %>%
    summarise(s_t_sq = var(!!sym(outcome)), .groups = "drop")

  s_j_sq_df <- matches_table %>%
    filter(!!sym(treatment) == 0) %>%
    left_join(s_t_sq_df, by = "subclass") %>%
    group_by(id) %>%
    summarise(s_j_sq = mean(s_t_sq, na.rm = TRUE), .groups = "drop")

  list(s_t_sq = s_t_sq_df, s_j_sq = s_j_sq_df)
}


#' Compute S1^2 via treated-to-treated matching
#'
#' For each treated unit t, finds the K nearest treated units, then
#' computes s_{1t}^2 = (Y_t - mean_{K-NN} Y)^2. Returns S1^2 =
#' mean_t(s_{1t}^2).
#'
#' @param df Full data frame (treated and control units).
#' @param treatment Name of treatment variable.
#' @param outcome Name of outcome variable.
#' @param K Number of nearest treated neighbors.
#' @param covs Character vector of covariate column names, or NULL to
#'   auto-detect columns starting with "X" (default NULL).
#' @param scaling Scaling vector passed to get_cal_matches.  If NULL,
#'   uses default_scaling(df, covs).
#' @param metric Distance metric passed to get_cal_matches (default
#'   "maximum").
#' @param id_name Name of the unit ID column (default "id").
#' @return A list with:
#'   \describe{
#'     \item{S1_sq}{Scalar: mean of s_{1t}^2 over treated units.}
#'     \item{s_1t_sq}{Data frame with columns \code{id} and \code{s_1t_sq}.}
#'   }
#' @export
calculate_S1_sq_treated_to_treated <- function(
    df, treatment = "Z", outcome = "Y",
    K = 1, covs = NULL,
    scaling = NULL, metric = "maximum", id_name = "id") {

  treatment_sym <- sym(treatment)
  id_sym <- sym(id_name)

  df <- df %>% mutate(!!id_sym := as.character(!!id_sym))
  df_t <- df %>% filter(!!treatment_sym == 1)

  # Resolve covariates to a character vector
  if (is.null(covs)) {
    covs <- get_x_vars(df)
  }

  df_t_pool <- df_t %>%
    mutate(!!treatment_sym := 0L,
           !!id_sym := paste0(!!id_sym, "_pool"))

  df_fake <- dplyr::bind_rows(df_t, df_t_pool)

  if (is.null(scaling)) {
    scaling <- default_scaling(df, covs)
  }

  matches_tt <- get_cal_matches(
    data = df_fake,
    covs = covs,
    treatment = treatment,
    metric = metric,
    rad_method = "knn",
    k = K + 1,
    scaling = scaling,
    id_name = id_name,
    warn = FALSE,
    est_method = "average"
  )

  s_1t_sq_list <- purrr::map(matches_tt$matches, function(match_df) {
    t_id <- as.character(match_df$subclass[1])
    t_Y <- df_t %>%
      filter(!!id_sym == t_id) %>%
      pull(!!sym(outcome))

    neighbors <- match_df %>%
      filter(!!treatment_sym == 0, dist > 1e-8) %>%
      slice_head(n = K)

    if (nrow(neighbors) < K) {
      warning(paste("Treated unit", t_id, "has only", nrow(neighbors),
                    "treated neighbours; s_1t_sq set to NA."))
      return(tibble(!!id_sym := t_id, s_1t_sq = NA_real_))
    }

    Y_hat <- mean(neighbors[[outcome]])
    s_1t_sq <- (t_Y - Y_hat)^2 * K / (K + 1L)
    tibble(!!id_sym := t_id, s_1t_sq = s_1t_sq)
  })

  s_1t_sq_df <- bind_rows(s_1t_sq_list)
  S1_sq <- mean(s_1t_sq_df$s_1t_sq, na.rm = TRUE)

  list(S1_sq = S1_sq, s_1t_sq = s_1t_sq_df)
}


#' Finite-sample variance estimator
#'
#' Implements the plug-in estimator
#'   hat_V = S1^2 / n_T  +  S0^2 / ESS(C)
#' where S0^2 is the w_j^2-weighted average of s_j^2 across control units, and
#' S1^2 is either the simple average of s_t^2 (common-variance assumption) or
#' the average of s_{1t}^2 from treated-to-treated K-NN matching.
#'
#' Also returns the estimated empirical covariance Cov_p(w_j, s_j^2)
#' which goes with the alternate formulation of:
#'
#'   hat_V_alt = S^2*(1/n_T + 1/ESS_C) + (1/n_T)*Cov_p(w_j, s_j^2)
#'
#' where S is a pooled estimate of S1 and S0
#'
#' @param matches Fitted CSM object (csm_matches) or a matched data frame
#'   (e.g. output of \code{full_unit_table()}).
#' @param df Full original data frame.  Required when
#'   \code{use_common_variance = FALSE} and \code{matches} is a data frame.
#' @param outcome Name of outcome variable (default "Y").
#' @param treatment Name of treatment variable (default "Z").
#' @param use_common_variance If TRUE (default), estimate S1^2 from the
#'   control-side subclass variances (assuming sigma_1(x)=sigma_0(x)).
#'   If FALSE, estimate S1^2 via treated-to-treated K-NN matching.
#' @param K Number of treated neighbours for the treated-to-treated step
#'   (used only when \code{use_common_variance = FALSE}).
#' @param covs Character vector of covariate names.  Ignored when
#'   \code{matches} is a CSM object (extracted from \code{params(matches)}).
#'   Required when \code{matches} is a data frame and
#'   \code{use_common_variance = FALSE}.
#' @param scaling Per-covariate scaling vector.  Same rules as \code{covs}.
#' @param metric Distance metric (default "maximum").  Same rules as \code{covs}.
#' @param id_name Name of the unit ID column (default "id").  Same rules as
#'   \code{covs}.
#'
#' @return A tibble with columns matching \code{get_total_variance()}:
#'   \code{V} (= V_E * N_T, for SE = sqrt(V)/sqrt(N_T) consistency),
#'   \code{V_E} (the finite-sample variance estimate S1^2/N_T + S0^2/ESS_C),
#'   \code{V_P} (NA — not decomposed by this estimator),
#'   \code{SE} (= sqrt(V_E)),
#'   \code{N_T}, \code{ESS_C}, \code{sigma_hat} (NA),
#'   plus diagnostics \code{S0_sq}, \code{S1_sq}, \code{cov_w_s}.
#' @export
get_finite_variance <- function(
    matches,
    df = NULL,
    outcome = "Y",
    treatment = "Z",
    use_common_variance = TRUE,
    homoskedastic = FALSE,
    K = 1,
    covs = NULL,
    scaling = NULL,
    metric = "maximum",
    id_name = "id") {

  # Unpack CSM object or accept a plain matched data frame
  if (is.csm_matches(matches)) {
    matches_table <- full_unit_table(matches) %>% mutate(id = as.character(id))
  } else {
    matches_table <- matches %>% mutate(id = as.character(id))
  }

  sample_sizes <- calc_N_T_N_C(matches_table, treatment = treatment)
  N_T   <- sample_sizes$N_T
  N_C <- sample_sizes$N_C
  ESS_C <- sample_sizes$N_C_tilde

  sjs <- calculate_s_j_sq(matches_table, outcome = outcome, treatment = treatment)
  p_drop = (N_T - nrow(sjs$s_t_sq)) / N_T
  s_t_sq_df <- sjs$s_t_sq
  s_j_sq_df <- sjs$s_j_sq

  if ( homoskedastic ) {
    V_E = mean( s_t_sq_df$s_t_sq ) * (1/N_T + 1/ESS_C)
    rs <- tibble(
      V_E       = V_E,
      SE        = sqrt(V_E),
      N_T       = N_T,
      N_C      = N_C,
      ESS_C     = ESS_C,
      sigma_hat = NA_real_,
      p_drop   = p_drop,
      S0_sq     = mean( s_t_sq_df$s_t_sq ),
      S1_sq    = S0_sq # under homoskedasticity, S1^2 = S0^2
    )
    return( rs )
  }

  w_j_df <- matches_table %>%
    filter(!!sym(treatment) == 0) %>%
    group_by(id) %>%
    summarise(w_j = sum(weights), .groups = "drop")

  control_df <- w_j_df %>%
    left_join(s_j_sq_df, by = "id")

  sum_w_sq <- sum(control_df$w_j^2, na.rm = TRUE)
  S0_sq    <- sum(control_df$w_j^2 * control_df$s_j_sq, na.rm = TRUE) / sum_w_sq
  if ( all( is.na( control_df$s_j_sq ) ) ) {
    cov_w_s = NA
  } else {
    cov_w_s <- cov(control_df$w_j, control_df$s_j_sq, use = "complete.obs")
  }

  if (use_common_variance) {
    S1_sq <- mean(s_t_sq_df$s_t_sq, na.rm = TRUE)
  } else {
    if ( is.csm_matches(matches) ) {
      df = full_unit_table(matches)
      pp          <- params(matches)
      covs    <- pp$covs
      scaling <- pp$scaling
      metric  <- pp$metric
      id_name <- pp$id_name
    }
    if (is.null(df)) {
      stop("'df' is required when use_common_variance = FALSE")
    }
    tt_result  <- calculate_S1_sq_treated_to_treated(
      df = df, treatment = treatment, outcome = outcome,
      K = K, covs = covs, scaling = scaling,
      metric = metric, id_name = id_name
    )
    S1_sq      <- tt_result$S1_sq
  }

  V_E <- S1_sq / N_T + S0_sq / ESS_C

  tibble(
    SE        = sqrt(V_E),
    N_T       = N_T,
    N_C       = N_C,
    ESS_C     = ESS_C,
    sigma_hat = NA_real_,
    V_E       = V_E,
    p_drop  = p_drop,
    S0_sq     = S0_sq,
    S1_sq     = S1_sq,
    cov_w_s   = cov_w_s
  )
}
