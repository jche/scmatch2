source("./scripts/inference-scripts/0_sim_inference_utils.R")
if ( FALSE ) {
  # --- Test run of one_iteration ---
  message("--- Testing one_iteration ---")
  needed_libs_iter <- c("tidyverse", "mvtnorm") # Add CSM if needed globally
  if (load_libs(needed_libs_iter)) {
    # Assuming necessary functions (gen_one_toy, etc.) are loaded
    # You might need to explicitly load CSM or its functions if not attached
    # e.g. library(CSM) or source("path/to/CSM_functions.R")

    source("./scripts/inference-scripts/0_sim_inference_utils.R")
    message("Running one_iteration with k=4...")
    iteration_result_k4 <- one_iteration( i = 1, k = 4, nc = 350, nt = 30, R = 1, verbose=TRUE )
    print(iteration_result_k4)

    iteration_result_k4 <- one_iteration( i = 1, k = 4, nc = 500, scaling = 8, toy_ctr_dist = 0.5, R = 1, verbose=TRUE )
    print(iteration_result_k4)

    message("\nRunning one_iteration with k=2...")
    source("./scripts/inference-scripts/0_sim_inference_utils.R")
    iteration_result_k2 <- one_iteration( i = 1, k = 2, nc = 100, nt = 25, R = 1, verbose=TRUE )
    print(iteration_result_k2)
  } else {
    message("Skipping one_iteration test due to missing packages.")
  }
}

if ( FALSE ) {
  # --- Test the save_match feature ---
  message("--- Testing save_match ---")
  needed_libs_iter <- c("tidyverse", "mvtnorm") # Add CSM if needed globally
  if (load_libs(needed_libs_iter)) {
    source("./scripts/inference-scripts/0_sim_inference_utils.R")
    iteration_result_k4 <- one_iteration( i = 1, k = 4, nc = 350, nt = 30, R = 1,
                                          save_match = TRUE, path_to_save_match = "./data/outputs/4d_bias_inference_Apr2025/test",
                                          verbose=TRUE )
  } else {
    message("Skipping one_iteration test due to missing packages.")
  }
}







# if ( FALSE ) {
#   # --- Example of running the full simulation for k=6 ---
#   message("\n--- Running full simulation example (k=6) ---")
#   needed_libs_sim <- c("tidyverse", "mvtnorm", "purrr") # Add CSM, furrr, future if needed
#   if (load_libs(needed_libs_sim)) {
#     # Set parallel = TRUE only if furrr and future are loaded and desired
#     use_parallel <- F
#     if (use_parallel) {
#       if (!requireNamespace("furrr", quietly = TRUE) || !requireNamespace("future", quietly = TRUE)) {
#         message("Packages 'furrr' and 'future' needed for parallel=TRUE. Running sequentially.")
#         use_parallel <- FALSE
#       } else {
#         library(furrr)
#         library(future)
#       }
#     }
#
#     # Ensure source files for data generation, matching, inference are loaded
#     # source("path/to/gen_toy_covar_k.R") etc.
#     source("./scripts/inference-scripts/0_sim_inference_utils.R")
#     globals <- list(
#       gen_one_toy = gen_one_toy,
#       gen_toy_covar_k = gen_toy_covar_k,
#       gen_df_adv_k = gen_df_adv_k,
#       `%||%` = `%||%`
#       # Add other custom functions used by gen_one_toy
#     )
#
#     results_k4 <- sim_inference_CSM_A_E( R = 2, # Smaller R for testing
#                                          k = 4,
#                                          nc = 500, # Smaller N for faster testing
#                                          nt = 75,
#                                          toy_ctr_dist = 0.5,
#                                          prop_nc_unif = 1/3,
#                                          # scaling = 10, # Adjust as needed
#                                          true_sigma = 0.5,
#                                          seed = 456,
#                                          parallel = use_parallel )
#     print(head(results_k4))
#     print(summary(results_k4))
#
#     # Example analysis (check coverage, bias)
#     if (requireNamespace("ggplot2", quietly = TRUE) && nrow(results_k4)>0) {
#       library(ggplot2)
#       results_k4 %>%
#         mutate(covered = (att_true >= CI_lower) & (att_true <= CI_upper)) %>%
#         summarize(k = first(k),
#                   R = n(),
#                   coverage = mean(covered, na.rm = TRUE),
#                   avg_bias_est = mean(att_est - att_true, na.rm = TRUE),
#                   avg_est_att = mean(att_est, na.rm=TRUE),
#                   avg_true_att = mean(att_true, na.rm=TRUE),
#                   avg_calc_bias = mean(bias, na.rm=TRUE), # Avg estimated bias post-matching
#                   avg_se = mean(se_AE, na.rm = TRUE),
#                   sd_att = sd(att_est, na.rm = TRUE), # Monte Carlo SD
#                   rmse = sqrt(mean( (att_est - att_true)^2, na.rm=TRUE))
#         ) -> sim_summary
#
#       print(sim_summary)
#
#       # Plot estimates vs true value
#       p1 <- ggplot(results_k4, aes(x = att_true, y = att_est)) +
#         geom_point(alpha = 0.5) +
#         geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
#         labs(title = paste("ATT Estimates vs True ATT (k=", unique(results_k4$k), ")"),
#              x = "True ATT", y = "Estimated ATT") +
#         theme_bw()
#       print(p1)
#
#       # Plot distribution of estimates
#       p2 <- ggplot(results_k4, aes(x=att_est)) +
#         geom_histogram(aes(y=after_stat(density)), bins=15, fill="lightblue", color="black", alpha=0.7) +
#         geom_vline(aes(xintercept=mean(att_true, na.rm=TRUE)), color="red", linetype="dashed", size=1) +
#         geom_density(color="blue", size=1) +
#         labs(title=paste("Distribution of ATT Estimates (k=", unique(results_k4$k), ")"),
#              x="Estimated ATT", y="Density",
#              caption = paste("Red line = Avg True ATT =", round(mean(results_k4$att_true, na.rm=TRUE), 3))) +
#         theme_bw()
#       print(p2)
#
#     }
#   } else {
#     message("Skipping full simulation example due to missing packages.")
#   }
# }
#
#
# if ( FALSE ) {
#   # --- Original diagnostic code block, updated ---
#   message("\n--- Running Diagnostic Code Block ---")
#   needed_libs_diag <- c("tidyverse", "mvtnorm", "ggplot2") # Add CSM if needed
#   if (load_libs(needed_libs_diag)) {
#
#     # Ensure source files for data generation, matching, inference are loaded
#     # library(CSM) or source("path/to/CSM_functions.R") etc.
#
#     k_dim <- 4 # Set dimension for diagnostics
#     nc_diag <- 500
#     nt_diag <- 125
#     message(paste("Generating diagnostic data (k=", k_dim, ")...", sep=""))
#
#     # Generate denoised data
#     df_dgp <- gen_one_toy( k = k_dim,
#                            nc = nc_diag,
#                            nt = nt_diag,
#                            ctr_dist=0.5,
#                            prop_nc_unif=1, # High overlap case
#                            f0_sd = 0) %>% # Generate denoised data
#       rename(Y0_denoised = Y0, Y1_denoised = Y1) %>%
#       mutate(Y_denoised = ifelse(Z, Y1_denoised, Y0_denoised)) %>%
#       select(-any_of("noise"))
#
#
#     # Summary of one dimension (e.g., X2 if k>=2)
#     if (k_dim >= 2) {
#       message("\nSummary of X2:")
#       print(summary( df_dgp$X2 ))
#     }
#
#     # Plot Y0 vs X1 colored by Treatment
#     p_diag1 <- ggplot( df_dgp, aes( X1, Y0_denoised, col=as.factor(Z) ) ) +
#       geom_point( alpha=0.2 ) +
#       geom_smooth( se=FALSE, method="loess", span=0.8 ) +
#       labs(title=paste("Y0 (denoised) vs X1 by Treatment (k=", k_dim,")", sep=""),
#            color="Treatment (Z)") +
#       theme_bw()
#     print(p_diag1)
#
#     # Check true treatment effect summary (denoised)
#     message("\nSummary of True Tx Effect (Y1_denoised - Y0_denoised):")
#     print(summary( df_dgp$Y1_denoised - df_dgp$Y0_denoised ))
#
#
#     scaling_diag <- 50
#     true_sigma_diag <- 0.5
#
#     # Add in the noise to set variation
#     df_dgp_i <- df_dgp %>%
#       mutate(noise = rnorm(n(), mean=0, sd=true_sigma_diag)) %>%
#       mutate(Y0 = Y0_denoised + noise,
#              Y1 = Y1_denoised + noise,
#              Y = ifelse(Z, Y1, Y0))
#     message("\nSummary of noisy data (df_dgp_i):")
#     print(head(df_dgp_i))
#
#
#     # Plot first two dimensions colored by treatment (if k >= 2)
#     if (k_dim >= 2) {
#       p_diag2 <- ggplot( df_dgp_i, aes( X1, X2, col=as.factor(Z) ) ) +
#         geom_point(alpha=0.5) +
#         labs(title=paste("X1 vs X2 by Treatment (k=",k_dim,")", sep=""),
#              color="Treatment (Z)") +
#         theme_bw()
#       print(p_diag2)
#     }
#
#     message("\nPerforming matching (k=", k_dim, ")...", sep="")
#     ### Perform matching
#     # CRITICAL ASSUMPTION: get_cal_matches handles k covariates correctly
#     mtch <- tryCatch(get_cal_matches(
#       df_dgp_i,
#       metric = "maximum",
#       scaling = scaling_diag,
#       caliper = 1,
#       rad_method = "adaptive-5nn",
#       est_method = "scm"
#       # Add covariate specification here if needed:
#       # , covariates = paste0("X", 1:k_dim)
#     ), error = function(e) { message("Matching failed: ", e$message); NULL})
#
#     if (!is.null(mtch)) {
#       print(mtch) # Print matching summary
#
#       # Can add plots using result_table or full_unit_table if needed and available
#
#       message("\nPerforming inference...")
#       ### Perform inference using the A-E method
#       ATT_estimate <- tryCatch(get_ATT_estimate( mtch ),
#                                error = function(e) { message("Inference failed: ", e$message); NULL})
#       print(ATT_estimate)
#
#       # Calculate results tibble for this single run (using one_iteration logic)
#       # Needs the full one_iteration code or refactoring to avoid repetition
#       # message("\nCalculating results for this run...")
#       # rs_diag <- one_iteration(i=1, k=k_dim, nc=nc_diag, nt=nt_diag, ...) # Use params from above
#       # print(rs_diag)
#
#     } else {
#       message("Skipping inference and further analysis as matching failed.")
#     }
#
#   } else {
#     message("Skipping diagnostic block due to missing packages.")
#   }
# }
#
