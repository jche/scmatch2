
# Utility functions called by "table_Section5.4.R" to calculate monte
# carlo standard errors.


kurtosis <- function(x){
  S_T = sd(x)
  kurt = mean( (x - mean(x))^4 ) / S_T^4
  return(kurt)
}

MCvar_SE <- function(x){
  S_T = sd(x); R <- length(x); k_T <- kurtosis(x)
  return(S_T^2 * sqrt( (k_T-1)/R ))
}

MCSE_SE <- function(x){
  return(sqrt(MCvar_SE(x)))
}

MCSE_bias <- function(x){
  S_T = sd(x); R <- length(x)
  return( sqrt(S_T^2 /R ))
}


#' Save results to file, possibly appending them if file exists.
save_res_to_csv<- function(curr_res,
                           FNAME) {

  if (file.exists(FNAME)) {
    write_csv(curr_res, FNAME, append=TRUE)
  } else {
    write_csv(curr_res, FNAME)
  }

} # save_res_to_csv




one_iteration <- function( i,
                           toy_ctr_dist=0.5,
                           scaling = 8,
                           true_sigma = 0.5,
                           prop_nc_unif = 1/3,
                           nc = 500,
                           verbose = TRUE ) {
  # i <- 1
  if( verbose ) {
    if ( i %% 10 == 0 || (i == R) ) {
      cat("iteration", i, "\n" )
    }
  }

  ### Make the data
  df_dgp <-
    gen_one_toy( nc = nc,
                 ctr_dist=toy_ctr_dist,
                 prop_nc_unif=prop_nc_unif) %>%
    mutate(Y0_denoised = Y0 - noise,
           Y1_denoised = Y1 - noise,
           Y_denoised = Y - noise)

  ## Add in the noise to set variation
  df_dgp_i <- df_dgp %>%
    mutate(noise = rnorm(n(), mean=0, sd=true_sigma)) %>%
    mutate(Y0 = Y0_denoised + noise,
           Y1 = Y1_denoised + noise,
           Y = Y_denoised + noise)


  ### Perform matching
  mtch <- get_cal_matches(
    df_dgp_i,
    metric = "maximum",
    scaling = scaling,
    caliper = 1,
    # rad_method = "adaptive",
    rad_method = "adaptive-5nn",
    est_method = "scm"
  )

  ### Perform inference using the A-E method
  ATT_estimate <- get_ATT_estimate( mtch )


  rs = tibble(
    runID = i,
    att_est = ATT_estimate$ATT,
    se_AE = ATT_estimate$SE,
    CI_lower = att_est -1.96 *  se_AE,
    CI_upper = att_est + 1.96 * se_AE,
    N_T = ATT_estimate$N_T,
    N_C_tilde = ATT_estimate$N_C_tilde,
    true_SE = true_sigma * sqrt(1 / N_T + 1 / N_C_tilde),
    CI_lower_with_true_SE = att_est - 1.96 * true_SE,
    CI_upper_with_true_SE = att_est + 1.96 * true_SE
  )

  # Get true ATT
  rs$att_true <-
    df_dgp_i %>%
    filter(Z ==1) %>%
    summarize(att = mean(Y1-Y0)) %>%
    pull(att)

  # Calculate true post-match bias, given the covariates.
  full_units <- full_unit_table(mtch, nonzero_weight_only = TRUE )
  rs$bias <- full_units %>%
    group_by(Z) %>%
    summarize(mn = sum(Y0_denoised*weights) / sum(weights)) %>%
    summarize(bias = last(mn) - first(mn)) %>%
    pull(bias)

  rs

}


if ( FALSE ) {
  one_iteration( 1 )

  library( tidyverse )
  library( CSM )

  df_dgp <- gen_one_toy( nc = 2000,
                 ctr_dist=0.5,
                 prop_nc_unif=1) %>%
    mutate(Y0_denoised = Y0 - noise,
           Y1_denoised = Y1 - noise,
           Y_denoised = Y - noise)

  summary( df_dgp$X2 )

  dd <- df_dgp %>%
    mutate( X2r = round( 4 * X2 ) )
  ggplot( dd, aes( X1, Y0_denoised, col=as.factor(Z) ) ) +
    facet_wrap( ~ X2r ) +
    geom_point( alpha=0.2 ) +
    geom_smooth( se=FALSE ) +
    theme_minimal()

  summary( df_dgp$Y1_denoised - df_dgp$Y0_denoised )

  scaling <- 50
  true_sigma <- 0.5

  ## Add in the noise to set variation
  df_dgp_i <- df_dgp %>%
    mutate(noise = rnorm(n(), mean=0, sd=true_sigma)) %>%
    mutate(Y0 = Y0_denoised + noise,
           Y1 = Y1_denoised + noise,
           Y = Y_denoised + noise)
  df_dgp_i

  ggplot( df_dgp_i, aes( X1, X2, col=as.factor(Z) ) ) +
    geom_point()


  ### Perform matching
  mtch <- get_cal_matches(
    df_dgp_i,
    metric = "maximum",
    scaling = 50,
    caliper = 1,
    # rad_method = "adaptive",
    rad_method = "adaptive-5nn",
    est_method = "scm"
  )
  mtch

  res <- result_table( mtch, return = "agg_co_units", feasible_only = FALSE )
  res
  ggplot( res, aes( X1, X2, color=as.factor(Z), size=weights ) ) +
    geom_point( alpha=0.5 ) +
    theme_minimal()

  full_units <- full_unit_table(mtch, nonzero_weight_only = FALSE )
  full_units
  sum( duplicated( full_units$id ) )

  ggplot( full_units, aes( X1, X2, color=as.factor(Z) ) ) +
    geom_point()

  ### Perform inference using the A-E method
  ATT_estimate <- get_ATT_estimate( mtch )
  ATT_estimate

  rs = tibble(
    runID = i,
    att_est = ATT_estimate$ATT,
    se_AE = ATT_estimate$SE,
    CI_lower = att_est -1.96 *  se_AE,
    CI_upper = att_est + 1.96 * se_AE,
    N_T = ATT_estimate$N_T,
    N_C_tilde = ATT_estimate$N_C_tilde,
    true_SE = true_sigma * sqrt(1 / N_T + 1 / N_C_tilde),
    CI_lower_with_true_SE = att_est - 1.96 * true_SE,
    CI_upper_with_true_SE = att_est + 1.96 * true_SE
  )

  # Get true ATT
  rs$att_true <-
    df_dgp_i %>%
    filter(Z ==1) %>%
    summarize(att = mean(Y1-Y0)) %>%
    pull(att)

  # Calculate true post-match bias, given the covariates.
  rs$bias <- full_units %>%
    group_by(Z) %>%
    summarize(mn = sum(Y0_denoised*weights) / sum(weights)) %>%
    summarize(bias = last(mn) - first(mn)) %>%
    pull(bias)

  rs
}


#' Run a simulation to check inference on the
#' A-E inference method
#'
#' @param dgp_name Name of the DGP
#' @param att0 True ATT
#' @param R Number of MC runs
#' @param toy_ctr_dist distance between centers in the to
#'
#' @return
#' @export
#'
#' @examples
sim_inference_CSM_A_E <- function(R=100,
                                  toy_ctr_dist=0.5,
                                  prop_nc_unif = 1/3,
                                  scaling = 8,
                                  true_sigma = 0.5,
                                  nc = 500,
                                  seed = NULL, parallel = FALSE ) {
  ### Example run inputs
  # R <- 10; toy_ctr_dist=0.5; dgp_name <- "toy"; att0<-F

  if ( is.null(seed) ) {
    warning( "Setting seed within sim_inference_CSM_A_E", call. = FALSE )
    set.seed(123)
  }

  if ( parallel ) {
    library( furrr )
    library( future )
    plan(multisession, workers = parallel::detectCores() - 1 )
    results = future_map_dfr( 1:R,
                              .f = one_iteration,
                              toy_ctr_dist = toy_ctr_dist,
                              prop_nc_unif = prop_nc_unif,
                              scaling = scaling,
                              true_sigma = true_sigma,
                              nc = nc,
                              verbose = FALSE,
                              .options = furrr_options(seed = NULL),
                              .progress = TRUE )
  } else {
    results <- map_df( 1:R, one_iteration,
                       toy_ctr_dist = toy_ctr_dist,
                       prop_nc_unif = prop_nc_unif,
                       scaling = scaling,
                       true_sigma = true_sigma,
                       nc = nc,
                       .progress = TRUE )
  }


  return(results)
}


# Testing
if ( FALSE ) {
  rs <- sim_inference_CSM_A_E( R = 10 )
  rs
}


if ( FALSE ) {
  # Some testing code to see overlap

  debug( sim_inference_CSM_A_E )
  sim_inference_CSM_A_E( R = 3, prop_nc_unif = 0.2, seed = 123, parallel = FALSE )

  sim_inference_CSM_A_E( R = 3, prop_nc_unif = 1, seed = 123, parallel = FALSE )
}


one_iteration_OR <- function( i,
                           toy_ctr_dist=0.5,
                           scaling = 8,
                           true_sigma = 0.5,
                           prop_nc_unif = 1/3,
                           nc = 500,
                           verbose = TRUE ) {
  # i <- 1
  if( verbose ) {
    if ( i %% 10 == 0 || (i == R) ) {
      cat("iteration", i, "\n" )
    }
  }

  ### Make the data
  df_dgp <-
    gen_one_toy( nc = nc,
                 ctr_dist=toy_ctr_dist,
                 prop_nc_unif=prop_nc_unif) %>%
    mutate(Y0_denoised = Y0 - noise,
           Y1_denoised = Y1 - noise,
           Y_denoised = Y - noise)

  ## Add in the noise to set variation
  df_dgp_i <- df_dgp %>%
    mutate(noise = rnorm(n(), mean=0, sd=true_sigma)) %>%
    mutate(Y0 = Y0_denoised + noise,
           Y1 = Y1_denoised + noise,
           Y = Y_denoised + noise)


  ### Perform matching
  mtch <- get_cal_matches(
    df_dgp_i,
    metric = "maximum",
    scaling = scaling,
    caliper = 1,
    # rad_method = "adaptive",
    rad_method = "adaptive-5nn",
    est_method = "scm"
  )

  # ### Perform inference using the A-E method
  # ATT_estimate <- get_ATT_estimate( mtch )

  ### Perform inference using the OR method
  get_ATT_estimate_OR <- function( scmatch, treatment = "Z", outcome = "Y" ) {

    ATT = get_att_ests( scmatch, treatment = treatment, outcome = outcome )
    se = get_se_OR( scmatch, treatment = treatment, outcome = outcome )
    # se = get_se_AE( scmatch, treatment = treatment, outcome = outcome )
    se$ATT = ATT

    se %>% relocate( ATT ) %>%
      mutate(  t = ATT/SE )
  }
  ATT_estimate <- get_ATT_estimate_OR( mtch )


  rs = tibble(
    runID = i,
    att_est = ATT_estimate$ATT,
    se_AE = ATT_estimate$SE,
    CI_lower = att_est -1.96 *  se_AE,
    CI_upper = att_est + 1.96 * se_AE,
    N_T = ATT_estimate$N_T,
    N_C_tilde = ATT_estimate$N_C_tilde,
    true_SE = true_sigma * sqrt(1 / N_T + 1 / N_C_tilde),
    CI_lower_with_true_SE = att_est - 1.96 * true_SE,
    CI_upper_with_true_SE = att_est + 1.96 * true_SE
  )

  # Get true ATT
  rs$att_true <-
    df_dgp_i %>%
    filter(Z ==1) %>%
    summarize(att = mean(Y1-Y0)) %>%
    pull(att)

  # Calculate true post-match bias, given the covariates.
  full_units <- full_unit_table(mtch, nonzero_weight_only = TRUE )
  rs$bias <- full_units %>%
    group_by(Z) %>%
    summarize(mn = sum(Y0_denoised*weights) / sum(weights)) %>%
    summarize(bias = last(mn) - first(mn)) %>%
    pull(bias)

  rs

}

# simulation inference using Otsu and Rai boostrap variance
#' Run a simulation to check inference on the
#' A-E inference method
#'
#' @param dgp_name Name of the DGP
#' @param att0 True ATT
#' @param R Number of MC runs
#' @param toy_ctr_dist distance between centers in the to
#'
#' @return
#' @export
#'
#' @examples
sim_inference_CSM_OR <- function(R=100,
                                  toy_ctr_dist=0.5,
                                  prop_nc_unif = 1/3,
                                  scaling = 8,
                                  true_sigma = 0.5,
                                  nc = 500,
                                  seed = NULL, parallel = FALSE ) {
  ### Example run inputs
  # R <- 10; toy_ctr_dist=0.5; dgp_name <- "toy"; att0<-F

  if ( is.null(seed) ) {
    warning( "Setting seed within sim_inference_CSM_A_E", call. = FALSE )
    set.seed(123)
  }

  if ( parallel ) {
    library( furrr )
    library( future )
    # plan(multisession, workers = parallel::detectCores() - 1 )
    plan(multisession, workers = min(4, parallel::detectCores() - 1))
    results = future_map_dfr( 1:R,
                              .f = one_iteration_OR,
                              toy_ctr_dist = toy_ctr_dist,
                              prop_nc_unif = prop_nc_unif,
                              scaling = scaling,
                              true_sigma = true_sigma,
                              nc = nc,
                              verbose = FALSE,
                              .options = furrr_options(seed = NULL),
                              .progress = TRUE )
  } else {
    results <- map_df( 1:R, one_iteration_OR,
                       toy_ctr_dist = toy_ctr_dist,
                       prop_nc_unif = prop_nc_unif,
                       scaling = scaling,
                       true_sigma = true_sigma,
                       nc = nc,
                       .progress = TRUE )
  }


  return(results)
}



# OLD VERSION ----


#' Run a simulation to check inference on the
#' A-E inference method
#'
#' @param dgp_name Name of the DGP
#' @param att0 True ATT
#' @param R Number of MC runs
#' @param toy_ctr_dist distance between centers in the to
#'
#' @return
#' @export
#'
#' @examples
sim_inference_CSM_A_E_OLD <- function(dgp_name,
                                      att0,
                                      R=100,
                                      toy_ctr_dist=0.5,
                                      prop_nc_unif = 1/3,
                                      seed = NULL) {
  ### Example run inputs
  # R <- 10; toy_ctr_dist=0.5; dgp_name <- "toy"; att0<-F
  covered <- CI_lower <- CI_upper <-
    covered_with_true_SE <- CI_lower_with_true_SE <- CI_upper_with_true_SE <-
    att_true <- att_est <- att_debiased <-
    se_AE <- true_SE <-
    error <- bias <-
    N_T <- N_C_tilde <- numeric(R)

  if ( is.null(seed) ) {
    set.seed(123)
    ### Generate one dataset
    df_dgp <-
      gen_one_toy(ctr_dist=toy_ctr_dist,
                  prop_nc_unif=prop_nc_unif) %>%
      mutate(Y0_denoised = Y0 - noise,
             Y1_denoised = Y1 - noise,
             Y_denoised = Y - noise)
  }
  scaling <- 8

  for (i in 1:R){
    # i <- 1
    print(i)

    if ( !is.null(seed) ) {
      set.seed(seed[i])
      ### Generate one dataset
      df_dgp <-
        gen_one_toy(ctr_dist=toy_ctr_dist,
                    prop_nc_unif=prop_nc_unif) %>%
        mutate(Y0_denoised = Y0 - noise,
               Y1_denoised = Y1 - noise,
               Y_denoised = Y - noise)
    }

    ## Re-generate the noises
    true_sigma <- 0.5
    df_dgp_i <- df_dgp %>%
      mutate(noise = rnorm(n(), mean=0, sd=true_sigma)) %>%
      mutate(Y0 = Y0_denoised + noise,
             Y1 = Y1_denoised + noise,
             Y = Y_denoised + noise)

    ### Perform matching
    mtch <- get_cal_matches(
      df_dgp_i,
      metric = "maximum",
      scaling = scaling,
      caliper = 1,
      # rad_method = "adaptive",
      rad_method = "adaptive-5nn",
      est_method = "scm"
    )

    ### Perform inference using the A-E method
    ATT_estimate <- get_ATT_estimate( mtch )
    att_est[i] <- ATT_estimate$ATT
    se_AE[i] = ATT_estimate$SE
    CI_lower[i] = att_est[i] -1.96 *  se_AE[i]
    CI_upper[i] = att_est[i] + 1.96 * se_AE[i]

    N_T[i] <- ATT_estimate$N_T
    N_C_tilde[i] <- ATT_estimate$N_C_tilde
    true_SE[i] = true_sigma *
      sqrt(1 / N_T[i] +1 / N_C_tilde[i])
    CI_lower_with_true_SE[i] =
      att_est[i] - 1.96 * true_SE[i]
    CI_upper_with_true_SE[i] =
      att_est[i] + 1.96 * true_SE[i]

    treatment_table <- mtch$treatment_table


    ### Obtain necessary things to output
    # Get true ATT
    if (att0){
      att_true[i] <- att <- 0
    }else{
      att_true[i] <- att <-
        df_dgp_i %>%
        filter(Z & (id %in% treatment_table$id)) %>%
        summarize(att = mean(Y1-Y0)) %>%
        pull(att)
    }

    covered[i] = (CI_lower[i] < att) & (att < CI_upper[i])
    covered_with_true_SE[i] =
      (CI_lower_with_true_SE[i] < att) &
      (att < CI_upper_with_true_SE[i])

    # Get error and bias
    full_units <- full_unit_table(mtch, nonzero_weight_only = TRUE )
    error[i] <- full_units %>%
      group_by(Z) %>%
      summarize(mn = sum(noise*weights) / sum(weights)) %>%
      summarize(est = last(mn) - first(mn)) %>%
      pull(est)

    bias[i] <- full_units %>%
      group_by(Z) %>%
      summarize(mn = sum(Y0_denoised*weights) / sum(weights)) %>%
      summarize(est = last(mn) - first(mn)) %>%
      pull(est)

    print(paste0("att_est is ", signif(att_est[i],3),
                 "; LB is ", signif(CI_lower[i],3),
                 "; UB is ", signif(CI_upper[i],3)))
    print(paste0("A-E s.e. is ", se_AE[i]))
    print(paste0("Covered is ", covered[i]))

  }

  res_save_bayesian_boot <-
    tibble(id=1:R,
           name=dgp_name,
           att_true = att_true,
           att_est= att_est,
           lower=CI_lower,
           upper=CI_upper,
           covered=covered,
           se_AE=se_AE,
           error=error,
           bias=bias,
           true_SE = true_SE,
           lower_with_true_SE=CI_lower_with_true_SE,
           upper_with_true_SE=CI_upper_with_true_SE,
           covered_with_true_SE=covered_with_true_SE,
           N_T = N_T,
           N_C_tilde = N_C_tilde)

  return(res_save_bayesian_boot)
}




# Testing ----

if ( FALSE ) {

  R = 3
  toy_naive_low <-
    sim_inference_CSM_A_E(
      dgp_name="toy",
      att0=F,
      R=R,
      prop_nc_unif = 1/3,
      seed = c(123 + 1:R * 2)
    )

}
