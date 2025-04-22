
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


#' Run one iteration of the simulation
#'
#' @param i Iteration number (integer).
#' @param k Dimensionality of covariates (integer, default 2).
#' @param toy_ctr_dist Control cluster distance for gen_one_toy (numeric).
#' @param scaling Scaling parameter for get_cal_matches (numeric).
#' @param true_sigma SD of the noise added to potential outcomes (numeric).
#' @param prop_nc_unif Proportion of uniform controls for gen_one_toy (numeric).
#' @param nc Number of control units (integer).
#' @param nt Number of treated units (integer, passed to gen_one_toy).
#' @param save_match Save the resulting match object? (logical, default FALSE). # <<< ADDED DOC
#' @param path_to_save_match Directory path to save the match object if save_match is TRUE (character). Filename will be generated based on iteration number. # <<< ADDED DOC
#' @param verbose Print progress message (logical).
#' @param R Total number of iterations (for verbose printing).
#'
#' @return A tibble with results for one iteration.
#' @importFrom stats rnorm weighted.mean reformulate
#' @importFrom dplyr %>% rename mutate select filter summarize pull group_by last first any_of n tibble
#' @importFrom rlang `%||%` # Assuming %||% comes from rlang or similar

one_iteration <- function( i,
                           k = 2,
                           toy_ctr_dist = 0.5,
                           scaling = 8,
                           true_sigma = 0.5,
                           prop_nc_unif = 1/3,
                           nc = 500,
                           nt = 100,
                           save_match = FALSE,         # <<< ADDED PARAM
                           path_to_save_match = "",    # <<< ADDED PARAM
                           verbose = TRUE,
                           R = NA ) {

  # Verbose output
  if( verbose ) {
    if (!is.na(R) && R > 0 && (i %% 10 == 0 || (i == R)) ) {
      cat("iteration", i, " (k=", k, ")\n", sep="")
    } else if (is.na(R) && i %% 10 == 0) { # Fallback if R not provided
      cat("iteration", i, " (k=", k, ")\n", sep="")
    }
  }

  ### Make the data using k dimensions
  df_dgp <-
    gen_one_toy( k = k,
                 nc = nc,
                 nt = nt,
                 ctr_dist=toy_ctr_dist,
                 prop_nc_unif=prop_nc_unif,
                 f0_sd = 0) %>%
    rename(Y0_denoised = Y0, Y1_denoised = Y1) %>%
    mutate(Y_denoised = ifelse(Z, Y1_denoised, Y0_denoised)) %>%
    select(-any_of("noise"))

  ## Add in the desired noise level for this iteration
  df_dgp_i <- df_dgp %>%
    mutate(noise = rnorm(n(), mean=0, sd=true_sigma)) %>%
    mutate(Y0 = Y0_denoised + noise,
           Y1 = Y1_denoised + noise,
           Y = ifelse(Z, Y1, Y0))

  ### Perform matching
  mtch <- tryCatch(
    get_cal_matches(
      df_dgp_i,
      metric = "maximum",
      scaling = scaling,
      caliper = 1,
      rad_method = "adaptive-5nn",
      est_method = "scm"
      # Consider adding covariate specification if needed:
      # , covariates = paste0("X", 1:k)
    ), error = function(e) {
      warning("Matching failed in iteration ", i, ": ", e$message, call. = FALSE)
      return(NULL)
    }
  )

  # --- >>> ADDED: Save the match object if requested <<< ---
  if (save_match && !is.null(mtch)) {
    # Check if the path is a non-empty string and the directory exists
    if (nzchar(path_to_save_match) && dir.exists(path_to_save_match)) {
      # Construct filename (e.g., match_iter_1.rds, match_iter_2.rds)
      file_name <- paste0(
        "match_iter_", i,
        "_k_", k,
        "_dist_", toy_ctr_dist,
        "_scale_", scaling,
        "_sigma_", true_sigma,
        "_prop_", round(prop_nc_unif, 4), # Round prop for cleaner filename
        "_nc_", nc,
        "_nt_", nt,
        ".rds"
      )
      full_path <- file.path(path_to_save_match, file_name)

      tryCatch({
        saveRDS(mtch, file = full_path)
        if (verbose) {
          cat("   -> Saved match object for iteration", i, "to:", full_path, "\n")
        }
      }, error = function(e) {
        # Warn if saving fails
        warning("Failed to save match object for iteration ", i, " to ", full_path, ": ", e$message, call. = FALSE)
      })
    } else if (nzchar(path_to_save_match)) {
      # Path provided but doesn't exist
      warning("Directory '", path_to_save_match, "' does not exist. Cannot save match object for iteration ", i, ".", call. = FALSE)
    } else {
      # Path is empty
      warning("'path_to_save_match' is empty. Cannot save match object for iteration ", i, ".", call. = FALSE)
    }
  } else if (save_match && is.null(mtch)) {
    # Optional: Notify if saving was requested but matching failed
    if (verbose) {
      cat("   -> Matching failed for iteration", i, ", skipping save.\n")
    }
  }
  # --- >>> END ADDED SECTION <<< ---


  ### Perform inference
  ATT_estimate <- NULL
  if (!is.null(mtch)) {
    ATT_estimate <- tryCatch(
      get_ATT_estimate( mtch ),
      error = function(e) {
        warning("ATT estimation failed in iteration ", i, ": ", e$message, call. = FALSE)
        return(NULL)
      }
    )
  }

  # Extract results safely
  att_est_val <- ATT_estimate$ATT %||% NA
  se_val <- ATT_estimate$SE %||% NA
  N_T_val <- ATT_estimate$N_T %||% NA
  ESS_C_val <- ATT_estimate$ESS_C %||% NA

  # Store iteration results
  rs = tibble(
    runID = i,
    k = k,
    att_est = att_est_val,
    SE = se_val,
    CI_lower = att_est_val - 1.96 * se_val, # Use val here
    CI_upper = att_est_val + 1.96 * se_val, # Use val here
    N_T = N_T_val,
    ESS_C = ESS_C_val
  )

  # Get true ATT (based on denoised potential outcomes for treated)
  rs$att_true <-
    df_dgp %>%
    filter(Z == TRUE) %>%
    summarize(att = mean(Y1_denoised - Y0_denoised)) %>%
    pull(att)

  # Calculate true post-match bias
  bias_val <- NA
  if (!is.null(mtch) && !is.null(ATT_estimate) && !is.na(att_est_val)) {
    # Try to get full units, including the necessary denoised outcome
    # NOTE: You might need to explicitly tell get_cal_matches or a subsequent
    # function to include Y0_denoised if it's not carried through automatically.
    # Assuming full_unit_table can access it or it's added manually.
    full_units <- full_unit_table(mtch, nonzero_weight_only = TRUE)

    # Proceed only if full_units is valid and has the required column
    if (!is.null(full_units) && nrow(full_units) > 0 && "Y0_denoised" %in% names(full_units) && "weights" %in% names(full_units)) {
      bias_calc <- full_units %>%
        group_by(Z) %>%
        summarize(mn = stats::weighted.mean(Y0_denoised, w = weights, na.rm=TRUE),
                  .groups = 'drop')

      if(nrow(bias_calc) == 2 && all(c(0, 1) %in% bias_calc$Z)) { # Check both groups exist
        # Bias = E[Y0|Z=1, M] - E[Y0|Z=0, M]
        mean_y0_t_matched <- bias_calc$mn[bias_calc$Z == 1]
        mean_y0_c_matched <- bias_calc$mn[bias_calc$Z == 0]
        bias_val <- mean_y0_t_matched - mean_y0_c_matched
      } else {
        warning("Bias calculation skipped in iteration ", i, ": Both Z groups not present or valid in weighted data.", call. = FALSE)
      }
    } else if (!is.null(full_units)) {
      if (!"Y0_denoised" %in% names(full_units)) warning("Bias calculation skipped in iteration ", i, ": Y0_denoised not found in full_units.", call. = FALSE)
      if (!"weights" %in% names(full_units)) warning("Bias calculation skipped in iteration ", i, ": 'weights' not found in full_units.", call. = FALSE)
    }
  } else {
    if (is.null(mtch)) warning("Bias calculation skipped in iteration ", i, ": Matching failed.", call. = FALSE)
    if (is.null(ATT_estimate)) warning("Bias calculation skipped in iteration ", i, ": ATT estimation failed.", call. = FALSE)
  }
  rs$bias <- bias_val

  return(rs)
}




#' Run a simulation study for CSM A-E inference
#'
#' @param R Number of Monte Carlo replications (integer).
#' @param k Dimensionality of covariates (integer, default 2).
#' @param toy_ctr_dist Control cluster distance for gen_one_toy (numeric).
#' @param prop_nc_unif Proportion of uniform controls for gen_one_toy (numeric).
#' @param scaling Scaling parameter for get_cal_matches (numeric).
#' @param true_sigma SD of the noise added to potential outcomes (numeric).
#' @param nc Number of control units (integer).
#' @param nt Number of treated units (integer).
#' @param seed Random seed to set for reproducibility (integer or NULL).
#' @param parallel Use parallel processing (logical or integer specifying workers).
#'
#' @return A tibble containing results from R iterations.
#' @export
#'
sim_inference_CSM_A_E <- function(R=100,
                                  k = 2, # <<< Added k parameter
                                  toy_ctr_dist=0.5,
                                  prop_nc_unif = 1/3,
                                  scaling = 8,
                                  true_sigma = 0.5,
                                  nc = 500,
                                  nt = 100, # <<< Make nt explicit
                                  seed = NULL,
                                  parallel = FALSE ) {

  # Input checks
  stopifnot(is.numeric(R), length(R) == 1, R >= 1)
  stopifnot(is.numeric(k), length(k) == 1, k >= 1)
  # ... add other checks as needed

  if ( !is.null(seed) ) {
    set.seed(seed)
    message(paste("Setting seed to", seed))
  } else {
    warning( "Seed not set for sim_inference_CSM_A_E run.", call. = FALSE )
  }

  # Prepare arguments list common to both parallel and sequential execution
  args_list <- list(
    k = k, # <<< Pass k
    toy_ctr_dist = toy_ctr_dist,
    prop_nc_unif = prop_nc_unif,
    scaling = scaling,
    true_sigma = true_sigma,
    nc = nc,
    nt = nt, # <<< Pass nt
    R = R # Pass R for verbose progress in one_iteration
  )

  run_start_time <- Sys.time()

  if ( identical(parallel, TRUE) || is.numeric(parallel) && parallel > 0) {
    # --- Parallel Execution ---
    if (!requireNamespace("furrr", quietly = TRUE) || !requireNamespace("future", quietly = TRUE)) {
      stop("Packages 'furrr' and 'future' are required for parallel execution.")
    }
    library(furrr)
    library(future)

    # Determine number of workers
    if (is.numeric(parallel)) {
      workers <- max(1, as.integer(parallel))
    } else {
      workers <- max(1, future::availableCores() - 1) # Default: N_cores - 1
    }
    plan(multisession, workers = workers)
    message(paste("Running", R, "iterations (k=", k,") in parallel on ", workers, " workers...", sep=""))


    # Use furrr_options(seed = TRUE) for automatic seedling based on main seed
    results = future_map_dfr( 1:R,
                              .f = one_iteration,
                              # Explicitly list each argument instead of using !!!args_list
                              k = k,
                              toy_ctr_dist = toy_ctr_dist,
                              prop_nc_unif = prop_nc_unif,
                              scaling = scaling,
                              true_sigma = true_sigma,
                              nc = nc,
                              nt = nt,
                              R = R,
                              # Other arguments
                              verbose = FALSE,
                              .options = furrr_options(
                                seed = TRUE,
                                scheduling = 1,
                                globals = globals
                              ),
                              .progress = TRUE )

    # Shut down workers after use
    plan(sequential)
    message("Parallel execution finished.")

  } else {
    # --- Sequential Execution ---
    if (!requireNamespace("purrr", quietly = TRUE)) {
      stop("Package 'purrr' is required for sequential execution with map_df.")
    }
    library(purrr)
    message(paste("Running", R, "iterations (k=", k, ") sequentially...", sep=""))

    # Check if progress package is available for the progress bar
    use_progress <- requireNamespace("progress", quietly = TRUE)

    results <- map_df( 1:R,
                       .f = one_iteration,
                       k = k,
                       toy_ctr_dist = toy_ctr_dist,
                       prop_nc_unif = prop_nc_unif,
                       scaling = scaling,
                       true_sigma = true_sigma,
                       nc = nc,
                       nt = nt,
                       R = R,
                       verbose = TRUE,
                       .progress = use_progress ) # Show progress bar if available
    message("Sequential execution finished.")
  }

  run_end_time <- Sys.time()
  message(paste("Total simulation time:", round(difftime(run_end_time, run_start_time, units="secs"), 2), "seconds"))

  return(results)
}


# --- Updated Testing / Example Blocks ---

# Utility function for safe requireNamespace/library loading
load_libs <- function(libs) {
  for (lib in libs) {
    if (!requireNamespace(lib, quietly = TRUE)) {
      warning(paste("Package", lib, "not found. Please install it."), call. = FALSE)
      return(FALSE)
    }
    library(lib, character.only = TRUE)
  }
  return(TRUE)
}

# `%||%` operator for safe default values
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (is.atomic(x) && all(is.na(x)))) y else x
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

    ATT = get_att_point_est( scmatch, treatment = treatment, outcome = outcome )
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
    ESS_C = ATT_estimate$ESS_C,
    true_SE = true_sigma * sqrt(1 / N_T + 1 / ESS_C),
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
    N_T <- ESS_C <- numeric(R)

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
    ESS_C[i] <- ATT_estimate$ESS_C
    true_SE[i] = true_sigma *
      sqrt(1 / N_T[i] +1 / ESS_C[i])
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
           ESS_C = ESS_C)

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
