#!/usr/bin/env Rscript

# --- Collection Script for 4D Bias Inference Simulation Results ---
# This script reads all result files from simulation iterations and combines them
# into a single stacked CSV file.

library(tidyverse)
library(here)

# --- Configuration ---
# Base directory (same as in the simulation script)
output_dir_base <- here::here("data/outputs/4d_bias_inference_Apr2025")
estimation_dir <- file.path(output_dir_base, "estimation")
combined_output_file <- file.path(output_dir_base, "combined_results.csv")

# --- Get list of all result files ---
cat("Looking for result files in:", estimation_dir, "\n")

result_files <- list.files(
  path = estimation_dir,
  pattern = "^results_iter_\\d+\\.rds$",
  full.names = TRUE
)

num_files <- length(result_files)

if (num_files == 0) {
  stop("No result files found in the directory. Check path and file pattern.", call. = FALSE)
} else {
  cat("Found", num_files, "result files.\n")
}

# --- Read and combine all results ---
cat("Reading and combining result files...\n")

# Progress tracking
total_files <- length(result_files)
progress_step <- max(1, floor(total_files / 10))  # Report progress at 10% increments

# Read all files and combine into a single dataframe
combined_results <- tibble()
failed_files <- character()

for (i in seq_along(result_files)) {
  file_path <- result_files[i]

  # Report progress
  if (i %% progress_step == 0 || i == total_files) {
    percent_done <- round(i / total_files * 100)
    cat(sprintf("  Progress: %d%% (%d/%d files)\n", percent_done, i, total_files))
  }

  # Try to read the file
  result <- tryCatch({
    readRDS(file_path)
  }, error = function(e) {
    file_name <- basename(file_path)
    warning("Failed to read file: ", file_name, " - ", e$message, call. = FALSE)
    failed_files <- c(failed_files, file_name)
    return(NULL)
  })

  # If successful, add to combined results
  if (!is.null(result)) {
    combined_results <- bind_rows(combined_results, result)
  }
}

# --- Report on the combined dataset ---
if (nrow(combined_results) == 0) {
  stop("No valid results were found. Check if files are properly formatted.", call. = FALSE)
} else {
  cat("\nSuccessfully combined", nrow(combined_results), "rows from",
      total_files - length(failed_files), "files.\n")

  if (length(failed_files) > 0) {
    cat("Failed to read", length(failed_files), "files.\n")
  }

  # Summary of iterations and overlap degrees
  iterations <- length(unique(combined_results$runID))
  overlap_degrees <- unique(combined_results$deg_overlap)

  cat("Dataset contains", iterations, "iterations across",
      length(overlap_degrees), "overlap degrees:",
      paste(overlap_degrees, collapse=", "), "\n")
}

# --- Save combined results to CSV ---
cat("Saving combined results to:", combined_output_file, "\n")

tryCatch({
  write_csv(combined_results, combined_output_file)
  cat("Combined results successfully saved.\n")
}, error = function(e) {
  stop("Failed to save the combined CSV file: ", e$message, call. = FALSE)
})

# --- Create a summary statistics file ---
summary_file <- file.path(output_dir_base, "results_summary.csv")
cat("Generating summary statistics...\n")

# Calculate summary statistics by overlap degree
summary_stats <- combined_results %>%
  mutate(ATT = 1.5 * k) %>%
  group_by(deg_overlap) %>%
  summarize(
    n_iterations = n_distinct(runID),
    mean_bias = mean(bias, na.rm = TRUE),
    sd_bias = sd(bias, na.rm = TRUE),
    mean_SATT = mean(SATT, na.rm = TRUE),
    mean_att_est = mean(att_est, na.rm = TRUE),
    mean_SE = mean(SE, na.rm = TRUE),
    coverage_rate = mean(CI_lower <= ATT & ATT <= CI_upper, na.rm = TRUE),
    coverage_rate_debiased = mean(CI_lower - bias <= ATT & ATT <= CI_upper-bias, na.rm = TRUE),
    mean_ESS_C = mean(ESS_C, na.rm = TRUE),
    mean_V_E = mean(V_E,na_rm=T),
    mean_sigma_hat = mean(sigma_hat, na.rm = TRUE),
    median_time_secs = median(time_secs, na.rm = TRUE),
    success_rate = mean(status == "Success", na.rm = TRUE)
  )

# Save summary statistics
tryCatch({
  write_csv(summary_stats, summary_file)
  cat("Summary statistics saved to:", summary_file, "\n")
}, error = function(e) {
  warning("Failed to save summary statistics: ", e$message, call. = FALSE)
})

cat("\nCollection process complete.\n")
