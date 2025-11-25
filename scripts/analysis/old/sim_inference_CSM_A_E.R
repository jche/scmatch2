


single_run <- function( toy_ctr_dist = 0.5,
                        prop_nc_unif = 1/3,
                        f0_sd = 0.5,
                        scaling = 8 ) {

  # toy_ctr_dist = 0.5; prop_nc_unif = 1/3, scaling = 8
  df_dgp <-
    gen_one_toy( ctr_dist=toy_ctr_dist,
                 prop_nc_unif=prop_nc_unif,
                 f0_sd = f0_sd )

  true_ATT = df_dgp %>%
    dplyr::filter( Z == 1 ) %>%
    summarise( ATT = mean( Y1 - Y0 ) ) %>%
    pull( ATT )

  ### Perform matching
  mtch <- get_cal_matches(
    df_dgp,
    metric = "maximum",
    scaling = scaling,
    caliper = 1,
    rad_method = "adaptive",
    est_method = "scm"
  )

  ### Perform inference using the A-E method
  ATT_estimate <- estimate_ATT( mtch ) %>%
    mutate( CI_lower = ATT - 1.96*SE,
            CI_upper = ATT + 1.96*SE,
            true_SE = f0_sd * sqrt( 1 / N_T + 1 / N_C_tilde ),
            CI_lower_true_SE = ATT - 1.96*true_SE,
            CI_upper_true_SE = ATT + 1.96*true_SE )



  ### Obtain necessary things to output

  # Get true ATT
  treatment_table <- mtch$treatment_table
  true_FATT <-
    df_dgp %>%
    filter( id %in% treatment_table$id ) %>%
    summarize(att = mean(Y1-Y0)) %>%
    pull(att)

  ATT_estimate <- mutate( ATT_estimate,
                          covered = (CI_lower < true_ATT) & (true_ATT < CI_upper),
                          covered_true_SE =
                            (CI_lower_true_SE < true_ATT) &
                            (true_ATT < CI_upper_true_SE) )

  # Get error and bias
  full_units <- result_table(mtch, nonzero_weight_only = TRUE )
  err_and_bias <- full_units %>%
    group_by(Z) %>%
    summarize(mn_noise = sum(noise*weights) / sum(weights),
              mn_f0 = sum((Y0-noise)*weights) / sum(weights) ) %>%
    summarize(err = last(mn_noise) - first(mn_noise),
              bias = last(mn_f0) - first(mn_f0) )

  ATT_estimate = bind_cols(ATT_estimate, err_and_bias)

  cat(paste0("att_est is ", signif(ATT_estimate$ATT,3),
               "; LB is ", signif(ATT_estimate$CI_lower,3),
               "; UB is ", signif(ATT_estimate$CI_upper,3)), "\n" )
  cat(paste0("\tA-E s.e. is ", signif( ATT_estimate$SE, 4 ) ), "\n" )
  cat(paste0("\tCovered is ", ATT_estimate$covered), "\n")

  ATT_estimate$trueATT = true_ATT
  ATT_estimate$trueFATT = true_FATT
  ATT_estimate
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
sim_inference_CSM_A_E <- function( dgp_name,
                                   att0,
                                   R=100,
                                   toy_ctr_dist=0.5,
                                   prop_nc_unif = 1/3,
                                   seed = NULL){
  ### Example run inputs
  # R <- 10; toy_ctr_dist=0.5; dgp_name <- "toy"; att0<-F

  if ( is.null(seed) ) {
    set.seed(123)
  }

  ### Generate one dataset
  reps = purrr::map_df( 1:R, function(r) {
    single_run( toy_ctr_dist = toy_ctr_dist,
                prop_nc_unif = prop_nc_unif,
                f0_sd = 0.5,
                scaling = 8 )
  } )

  return(reps)

}


summarize_sim_results <- function( reps ) {
  reps %>%
    summarize( coverage = mean(covered),
               coverage_true = mean(covered_true_SE),
               ESEhat = mean(SE),
               SE = sd( ATT ),
               bias = mean( ATT - trueATT ),
               RMSE = sqrt( mean( (ATT - trueATT)^2 ) ),
               E_NC = mean( N_C_tilde ) )
}

if ( FALSE ) {
  reps <- sim_inference_CSM_A_E( "toy", att0 = FALSE, R = 100 )
  summarize_sim_results( reps )
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
sim_inference_CSM_A_E_finite <- function( dgp_name,
                                          att0,
                                          R=100,
                                          toy_ctr_dist=0.5,
                                          prop_nc_unif = 1/3,
                                          seed = NULL){
  ### Example run inputs
  # R <- 10; toy_ctr_dist=0.5; dgp_name <- "toy"; att0<-F

  covered <- CI_lower <- CI_upper <-
    covered_true_SE <- CI_lower_true_SE <- CI_upper_true_SE <-
    att_true <- att_est <- att_debiased <-
    se_AE <- true_SE <-
    error <- bias <-
    N_T <- N_C_tilde <- numeric(R)

  if ( is.null(seed) ) {
    set.seed(123)
  }

  ### Generate one dataset
  df_dgp <-
    gen_one_toy(ctr_dist=toy_ctr_dist,
                prop_nc_unif=prop_nc_unif) %>%
    mutate(Y0_denoised = Y0 - noise,
           Y1_denoised = Y1 - noise,
           Y_denoised = Y - noise)
  scaling <- 8

  for (i in 1:R){
    # i <- 1
    print(i)

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
      rad_method = "adaptive",
      est_method = "scm"
    )

    ### Perform inference using the A-E method
    ATT_estimate <- estimate_ATT( mtch )
    att_est[i] <- ATT_estimate$ATT
    se_AE[i] = ATT_estimate$SE
    CI_lower[i] = att_est[i] -1.96 *  se_AE[i]
    CI_upper[i] = att_est[i] + 1.96 * se_AE[i]

    N_T[i] <- ATT_estimate$N_T
    N_C_tilde[i] <- ATT_estimate$N_C_tilde
    true_SE[i] = true_sigma *
      sqrt(1 / N_T[i] +1 / N_C_tilde[i])
    CI_lower_true_SE[i] =
      att_est[i] - 1.96 * true_SE[i]
    CI_upper_true_SE[i] =
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
    covered_true_SE[i] =
      (CI_lower_true_SE[i] < att) &
      (att < CI_upper_true_SE[i])

    # Get error and bias
    full_units <- result_table(mtch, nonzero_weight_only = TRUE )
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
           lower_true_SE=CI_lower_true_SE,
           upper_true_SE=CI_upper_true_SE,
           covered_true_SE=covered_true_SE,
           N_T = N_T,
           N_C_tilde = N_C_tilde)
  return(res_save_bayesian_boot)
}

