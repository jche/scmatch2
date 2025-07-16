source(here::here("scripts/inference-scripts/0_sim_inference_utils.R"))
# -------------------
# --- Simulation Parameters ---
# -------------------

# --- Main parameters: n, degree of freedom, error distribution ---
param_grid <- expand.grid(
  N = c(120, 240, 600, 1200),
  overlap = c("very_low", "low", "mid", "high", "very_high"),
  # error = c("homoskedastic", "covariate_dep", "treatment_dep", "location_dep"),
  error = c("homoskedastic", "covariate_dep", "treatment_dep"),
  stringsAsFactors = FALSE
)

# --- Supporting parameters: Fix them for a fixed k_dim. With differnt k_dim we might want to change ---
k_dim <- 4
toy_ctr_dist_base <- 0.5
scaling_base <- 8
include_bootstrap_global <- TRUE
boot_mtd_global <- "wild"
B_global <- 250

# --- Output Directories ---
paths <- get_sim_paths()
individual_dir <- paths$individual_dir
dir.create(individual_dir, showWarnings = FALSE, recursive = TRUE)

# --- Argument Input ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No iteration number supplied. Usage: Rscript parallel_sim_inference.R <iteration_number>", call. = FALSE)
} else {
  iteration_i <- as.integer(args[1])
  if (is.na(iteration_i)) {
    stop("Iteration number must be an integer.", call. = FALSE)
  }
}

# -------------------
# --- Main script ---
# -------------------
library(tidyverse)
library(mvtnorm)
library(doParallel)
library(foreach)
library(CSM)


cat("--- Starting Iteration:", iteration_i, "---\n")
set.seed(123 + iteration_i * 10) # Ensure reproducibility for each iteration


n_cores <- 60
cl <- makeCluster(n_cores)
registerDoParallel(cl)

results <- foreach(i = 1:nrow(param_grid), .combine = dplyr::bind_rows,
                   .packages = c("dplyr", "tibble", "purrr", "mvtnorm","CSM")) %dopar% {
  row <- param_grid[i, ]
  sim_master(
    iteration = iteration_i,
    N = row$N,
    overlap_label = row$overlap,
    error_label = row$error,
    grid_id = i,
    # Pass all supporting variables explicitly
    toy_ctr_dist_base = 0.5,
    scaling_base = 8,
    include_bootstrap_global = TRUE,
    boot_mtd_global = "wild",
    B_global = 250,
    k_dim = 4
  )
}


stopCluster(cl)

# --- Combine and Save Results for Iteration i ---
estimation_filename <- paste0("results_iter_", iteration_i, ".rds")
estimation_filepath <- file.path(individual_dir, estimation_filename)
saveRDS(results, file = estimation_filepath)
