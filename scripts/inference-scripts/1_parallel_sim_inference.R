library(tidyverse)
require(mvtnorm)
devtools::load_all()

source(here::here("scripts/inference-scripts/0_sim_inference_utils.R"))

# --- Simulation Parameters ---
# --- Main parameters: n, degree of freedom, error distribution ---
nc <- 500
nc_to_nt <- 5 / 1
nt <- nc / nc_to_nt

prop_nc_unif_values <- c(1/10, 1/5, 1/3, 1/2, 2/3)
deg_overlap_labels <- c("very_low", "low", "mid", "high", "very_high")

# --- Supporting parameters: Fix them ---
k_dim <- 4
toy_ctr_dist_base <- 0.5
scaling_base <- 8
true_sigma_base <- 0.5
include_bootstrap_global <- TRUE
boot_mtd_global <- "wild"
B_global <- 250

# --- Output Directories ---
# output_dir_base <- here::here("data/outputs/4d_bias_inference_Apr2025")
# individual_dir <- file.path(output_dir_base, "estimation")
paths <- get_sim_paths()
individual_dir <- paths$individual_dir

dir.create(individual_dir, showWarnings = FALSE, recursive = TRUE)

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

iteration_results_list <- list()

for (j in seq_along(prop_nc_unif_values)) {
  prop_unif <- prop_nc_unif_values[j]
  overlap_label <- deg_overlap_labels[j]

  cat("   Processing overlap:", overlap_label, "(prop_unif =", round(prop_unif, 3), ")\n")
  start_time <- Sys.time()

  # --- Call one_iteration Function ---
  single_overlap_result <- tryCatch({
    one_iteration(
      i = iteration_i,
      k = k_dim,
      nc = nc,
      nt = nt,
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
      nc = nc,               # Add nc here
      nt = nt,               # Add nt here
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
        nc = if (!"nc" %in% names(.)) nc else nc,
        nt = if (!"nt" %in% names(.)) nt else nt,
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
valid_results <- Filter(Negate(is.null), iteration_results_list)

if (length(valid_results) > 0) {
  final_results_tibble <- bind_rows(valid_results)

  estimation_filename <- paste0("results_iter_", iteration_i, ".rds")
  estimation_filepath <- file.path(individual_dir, estimation_filename)

  tryCatch({
    saveRDS(final_results_tibble, file = estimation_filepath)
    cat("--- Iteration", iteration_i, "Complete. Results saved to:", estimation_filepath, "---\n")
  }, error = function(e) {
    warning("Failed to save final results for iteration ", iteration_i, ": ", e$message, call. = FALSE)
  })
} else {
  cat("--- Iteration", iteration_i, "Failed: No valid results generated. ---\n")
}

# Clean exit
q(save = "no", status = 0)
