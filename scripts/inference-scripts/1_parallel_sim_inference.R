# Load necessary libraries
library(tidyverse)
require(mvtnorm)
# Assuming 'load_all()' loads the necessary custom functions like gen_one_toy, get_cal_matches, get_ATT_estimate
# If CSM is a package, use library(CSM) instead of load_all()
devtools::load_all() # Uncomment or replace with library(CSM) if needed

# Source the utility functions, including the assumed modified 'one_iteration'
source(here::here("scripts/inference-scripts/0_sim_inference_utils.R")) # Ensure this path is correct

# --- Simulation Parameters (Set Globally) ---
k_dim <- 4 # Covariate dimension as requested
nc_base <- 500 # Number of control units
nt_base <- 100 # Number of treated units
toy_ctr_dist_base <- 0.5 # Control cluster distance
scaling_base <- 8 # Scaling parameter for matching
true_sigma_base <- 0.5 # SD of noise added to potential outcomes
prop_nc_unif_values <- c(1/10, 1/5, 1/3, 1/2, 2/3) # 5 degrees of overlap
deg_overlap_labels <- c("very_low", "low", "mid", "high", "very_high") # Labels for overlap

include_bootstrap_global <- TRUE
boot_mtd_global <- "wild"
B_global <- 250

# --- Output Directories ---
# NOTE: one_iteration might handle saving internally, or might not save matches.
# This script no longer explicitly saves match objects. Adjust if needed.
output_dir_base <- here::here("data/outputs/4d_bias_inference_Apr2025")
# match_dir <- file.path(output_dir_base, "match") # Removed, assuming one_iteration handles or skips
estimation_dir <- file.path(output_dir_base, "estimation")

# Create directories if they don't exist
# dir.create(match_dir, showWarnings = FALSE, recursive = TRUE) # Removed
dir.create(estimation_dir, showWarnings = FALSE, recursive = TRUE)

# --- Get Iteration Number from Command Line ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No iteration number supplied. Usage: Rscript parallel_sim_inference.R <iteration_number>", call. = FALSE)
} else {
  iteration_i <- as.integer(args[1])
  if (is.na(iteration_i)) {
    stop("Iteration number must be an integer.", call. = FALSE)
  }
}

cat("--- Starting Iteration:", iteration_i, "---\n")
set.seed(123 + iteration_i * 10) # Ensure reproducibility for each iteration

# --- Store results for this iteration ---
iteration_results_list <- list()

# --- Loop Over Degrees of Overlap ---
for (j in seq_along(prop_nc_unif_values)) {
  prop_unif <- prop_nc_unif_values[j]
  overlap_label <- deg_overlap_labels[j]

  cat("   Processing overlap:", overlap_label, "(prop_unif =", round(prop_unif, 3), ")\n")
  start_time <- Sys.time()

  # --- Call one_iteration Function ---
  # Assumes one_iteration performs data generation, matching, inference,
  # and returns a tibble including columns: SATT, bias, att_est, SE, CI_lower, CI_upper, N_T, ESS_C etc.
  # It also needs to accept all the relevant parameters.
  # IMPORTANT: Ensure 'one_iteration' in '0_sim_inference_utils.R' is modified accordingly.
  single_overlap_result <- tryCatch({
    one_iteration(
      i = iteration_i,
      k = k_dim,
      nc = nc_base,
      nt = nt_base,
      toy_ctr_dist = toy_ctr_dist_base,
      prop_nc_unif = prop_unif,
      scaling = scaling_base,
      true_sigma = true_sigma_base,
      include_bootstrap = include_bootstrap_global,
      boot_mtd = boot_mtd_global,
      B = B_global,
      seed_addition = iteration_i,
      verbose = FALSE
    )
  }, error = function(e) {
    warning("one_iteration failed for iteration ", iteration_i, ", overlap ", overlap_label, ": ", e$message, call. = FALSE)
    # Return a placeholder tibble on error to maintain structure
    return(tibble(
      runID = iteration_i, # Assuming one_iteration would add this, but add here for consistency
      k = k_dim,
      deg_overlap = overlap_label, # Add overlap label here
      prop_nc_unif = prop_unif,   # Add prop_unif here
      nc = nc_base,               # Add nc here
      nt = nt_base,               # Add nt here
      toy_ctr_dist = toy_ctr_dist_base, # Add toy_ctr_dist here
      scaling = scaling_base,         # Add scaling here
      true_sigma = true_sigma_base,   # Add true_sigma here
      SATT = NA_real_,            # Placeholder for Sample ATT
      att_est = NA_real_,
      SE = NA_real_,
      CI_lower = NA_real_,
      CI_upper = NA_real_,
      bias = NA_real_,            # Placeholder for bias
      N_T = NA_integer_,
      ESS_C = NA_real_,
      V_E = NA_real_,
      sigma_hat = NA_real_,
      status = "Iteration Failed" # Add a status column
      # Add other columns returned by one_iteration as NA if needed
    ))
  }
  )

  # --- Record Time and Store Results ---
  end_time <- Sys.time()
  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat("      Overlap", overlap_label, "took", round(duration, 2), "seconds.\n")

  # Add timing information and potentially missing columns if one_iteration failed
  # Ensure runID, deg_overlap, prop_nc_unif etc. are present
  if (!is.null(single_overlap_result)) {
    single_overlap_result <- single_overlap_result %>%
      mutate(
        runID = iteration_i, # Ensure runID is set
        deg_overlap = overlap_label, # Ensure overlap label is set
        prop_nc_unif = prop_unif, # Ensure prop_unif is set
        time_secs = duration, # Add the calculated duration
        # Add other parameters if not returned by one_iteration
        k = if (!"k" %in% names(.)) k_dim else k,
        nc = if (!"nc" %in% names(.)) nc_base else nc,
        nt = if (!"nt" %in% names(.)) nt_base else nt,
        toy_ctr_dist = if (!"toy_ctr_dist" %in% names(.)) toy_ctr_dist_base else toy_ctr_dist,
        scaling = if (!"scaling" %in% names(.)) scaling_base else scaling,
        true_sigma = if (!"true_sigma" %in% names(.)) true_sigma_base else true_sigma,
        status = if (!"status" %in% names(.)) "Success" else status # Add status if missing
      ) %>%
      # Select and order columns for consistency
      select(
        runID, k, deg_overlap, prop_nc_unif, nc, nt, toy_ctr_dist, scaling, true_sigma,
        SATT, att_est, SE, CI_lower, CI_upper, bias, N_T, ESS_C, V_E, sigma_hat, # Assumed columns from one_iteration
        time_secs, status, everything() # Add time, status, and any other columns
      )
    iteration_results_list[[overlap_label]] <- single_overlap_result
  }

} # End loop over overlap degrees

# --- Combine and Save Results for Iteration i ---
# Filter out NULLs in case an error occurred before assignment
valid_results <- Filter(Negate(is.null), iteration_results_list)

if (length(valid_results) > 0) {
  final_results_tibble <- bind_rows(valid_results)

  estimation_filename <- paste0("results_iter_", iteration_i, ".rds")
  estimation_filepath <- file.path(estimation_dir, estimation_filename)

  tryCatch({
    saveRDS(final_results_tibble, file = estimation_filepath)
    cat("--- Iteration", iteration_i, "Complete. Results saved to:", estimation_filepath, "---\n")
  }, error = function(e) {
    warning("Failed to save final results for iteration ", iteration_i, ": ", e$message, call. = FALSE)
  })
} else {
  cat("--- Iteration", iteration_i, "Failed: No valid results generated. ---\n")
  # Optionally save an empty file or log the failure more formally
}

# Clean exit
q(save = "no", status = 0)
