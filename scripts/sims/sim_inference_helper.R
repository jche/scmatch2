save_res_to_csv<-
  function(curr_res,
           FNAME){
    if (file.exists(FNAME)) {
      write_csv(curr_res, FNAME, append=TRUE)
    } else {
      write_csv(curr_res, FNAME)
    }
  } # save_res_to_csv

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
sim_inference_CSM_A_E <- function(dgp_name,
                                  att0,
                                  R=100,
                                  toy_ctr_dist=0.5,
                                  prop_nc_unif = 1/3,
                                  seed = NULL){
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
    att_est[i] <- ATT_estimate$ATE
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

