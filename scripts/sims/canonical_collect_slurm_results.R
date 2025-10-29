# scripts/analysis/collect_slurm_results.R
# Collect and combine results from SLURM array jobs

library(tidyverse)

#' Collect results from a simulation
#'
#' @param sim_type One of "acic", "hainmueller", "kang"
#' @param n_expected Expected number of iterations
collect_sim_results <- function(sim_type, n_expected = 100) {

  cat(sprintf("\n=== Collecting %s results ===\n", toupper(sim_type)))

  # Directory with iteration results
  results_dir <- file.path("data", "outputs", "sim_slurm", sim_type)

  if (!dir.exists(results_dir)) {
    stop(sprintf("Results directory not found: %s", results_dir))
  }

  # Find all result files
  result_files <- list.files(
    results_dir,
    pattern = "^iter_\\d{4}\\.csv$",
    full.names = TRUE
  )

  cat(sprintf("Found %d result files (expected %d)\n",
              length(result_files), n_expected))

  if (length(result_files) == 0) {
    stop("No result files found!")
  }

  # Read and combine
  cat("Reading and combining results...\n")
  results <- map_df(result_files, read_csv, show_col_types = FALSE)

  # Check for duplicates
  if (any(duplicated(results$runid))) {
    warning("Found duplicate runids!")
  }

  # Check for missing iterations
  missing_ids <- setdiff(1:n_expected, results$runid)
  if (length(missing_ids) > 0) {
    warning(sprintf("Missing %d iterations: %s",
                    length(missing_ids),
                    paste(head(missing_ids, 10), collapse = ", ")))
  }

  # Save combined results
  output_file <- file.path("data", "outputs", "sim_slurm",
                           paste0(sim_type, "_combined.csv"))
  write_csv(results, output_file)

  cat(sprintf("Combined results saved to: %s\n", output_file))
  cat(sprintf("Total rows: %d\n", nrow(results)))
  cat(sprintf("Total columns: %d\n", ncol(results)))

  # Print summary
  cat("\nPerformance summary:\n")
  summary_stats <- results %>%
    pivot_longer(diff:kbal, names_to = "method", values_to = "estimate") %>%
    group_by(method) %>%
    summarize(
      rmse = sqrt(mean((estimate - true_ATT)^2, na.rm = TRUE)),
      bias = mean(estimate - true_ATT, na.rm = TRUE),
      n_na = sum(is.na(estimate)),
      .groups = "drop"
    ) %>%
    arrange(rmse)

  print(summary_stats, n = Inf)

  cat("\nTiming summary:\n")
  timing_stats <- results %>%
    summarize(
      mean_time = mean(elapsed_time_secs),
      median_time = median(elapsed_time_secs),
      min_time = min(elapsed_time_secs),
      max_time = max(elapsed_time_secs),
      total_time_hours = sum(elapsed_time_secs) / 3600
    )
  print(timing_stats)

  return(results)
}

# Collect all simulations
cat("\n", rep("=", 80), "\n", sep = "")
cat("COLLECTING ALL SIMULATION RESULTS\n")
cat(rep("=", 80), "\n", sep = "")

res_acic <- collect_sim_results("acic", n_expected = 100)
res_hain <- collect_sim_results("hainmueller", n_expected = 100)
res_kang <- collect_sim_results("kang", n_expected = 100)

cat("\n", rep("=", 80), "\n", sep = "")
cat("ALL RESULTS COLLECTED\n")
cat(rep("=", 80), "\n", sep = "")
