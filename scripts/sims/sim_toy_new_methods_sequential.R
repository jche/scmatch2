# Simulation script: scripts/sims/sim_toy_new_methods_sequential.R
# Run the three new methods sequentially with detailed timing

require(tidyverse)
require(mvtnorm)
require(grf)
require(twang)
require(kbal)
require(tictoc)
library(CSM)

# Create output directory
dir.create("data/outputs/sim_toy_results/toy_new_methods_seq",
           showWarnings = FALSE, recursive = TRUE)

source(here::here("R/wrappers.R"))
source(here::here("R/utils.R"))

# Number of iterations
n_iter <- 60

# Store all results
all_results <- list()

for (i in 1:n_iter) {

  cat("\n========================================\n")
  cat(sprintf("Starting iteration %d of %d\n", i, n_iter))
  cat("========================================\n")

  # Use same seed as main simulation
  set.seed(1000 + i)

  iter_start <- Sys.time()

  # ========================================================================
  # GENERATE DATA (same as main simulation)
  # ========================================================================
  cat("Generating data...\n")
  data_start <- Sys.time()

  nc <- 500
  nt <- 100
  f0_sd <- 0.5

  df <- gen_df_adv(
    nc = nc,
    nt = nt,
    f0_sd = f0_sd,
    tx_effect_fun = function(X1, X2) {3 * X1 + 3 * X2},
    f0_fun = function(x, y) {
      matrix(c(x, y), ncol = 2) %>%
        dmvnorm(
          mean = c(0.5, 0.5),
          sigma = matrix(c(1, 0.8, 0.8, 1), nrow = 2)
        ) * 20
    }
  )

  num_bins <- 6
  dist_scaling <- df %>%
    summarize(across(
      starts_with("X"),
      function(x) {
        if (is.numeric(x)) 2 * num_bins / (max(x) - min(x))
        else 1000
      }
    ))

  # Apply same filtering as main simulation
  preds_csm <- get_cal_matches(
    df = df,
    metric = "maximum",
    scaling = dist_scaling,
    rad_method = "fixed",
    est_method = "average",
    k = 25
  )

  df <- df %>%
    filter(!id %in% attr(preds_csm, "unmatched_units"))

  data_time <- as.numeric(difftime(Sys.time(), data_start, units = "secs"))
  cat(sprintf("  Data generation: %.2f seconds\n", data_time))

  my_covs <- c("X1", "X2")
  zform2 <- as.formula("Z ~ X1+X2")

  # ========================================================================
  # COMPUTE TRUE ATT
  # ========================================================================
  true_ATT <- df %>%
    filter(Z) %>%
    summarize(att = mean(Y1 - Y0)) %>%
    pull(att)

  cat(sprintf("  True ATT: %.3f\n", true_ATT))

  # ========================================================================
  # COMPUTE ESTIMATES - NEW METHODS
  # ========================================================================

  # Causal Forest
  cat("\nCausal Forest...\n")
  cf_start <- Sys.time()
  att_cf <- tryCatch(
    {
      result <- get_att_causal_forest(df, covs = my_covs)
      cat(sprintf("  Estimate: %.3f\n", result))
      result
    },
    error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
      return(NA_real_)
    }
  )
  cf_time <- as.numeric(difftime(Sys.time(), cf_start, units = "secs"))
  cat(sprintf("  Time: %.2f seconds\n", cf_time))

  # TWANG
  cat("\nTWANG...\n")
  twang_start <- Sys.time()
  att_twang <- tryCatch(
    {
      result <- get_att_twang(df, form = zform2)
      cat(sprintf("  Estimate: %.3f\n", result))
      result
    },
    error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
      return(NA_real_)
    }
  )
  twang_time <- as.numeric(difftime(Sys.time(), twang_start, units = "secs"))
  cat(sprintf("  Time: %.2f seconds\n", twang_time))

  # Kbal
  cat("\nKbal...\n")
  kbal_start <- Sys.time()
  att_kbal <- tryCatch(
    {
      result <- get_att_kbal(df, covs = my_covs, numdims = 30)
      cat(sprintf("  Estimate: %.3f\n", result))
      result
    },
    error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
      return(NA_real_)
    }
  )
  kbal_time <- as.numeric(difftime(Sys.time(), kbal_start, units = "secs"))
  cat(sprintf("  Time: %.2f seconds\n", kbal_time))

  # ========================================================================
  # CALCULATE TOTAL TIMING
  # ========================================================================
  total_time <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))

  cat("\n----------------------------------------\n")
  cat(sprintf("Iteration %d complete!\n", i))
  cat(sprintf("  Total time: %.2f seconds\n", total_time))
  cat(sprintf("  Breakdown:\n"))
  cat(sprintf("    Data gen:       %.2f sec (%.1f%%)\n",
              data_time, 100*data_time/total_time))
  cat(sprintf("    Causal Forest:  %.2f sec (%.1f%%)\n",
              cf_time, 100*cf_time/total_time))
  cat(sprintf("    TWANG:          %.2f sec (%.1f%%)\n",
              twang_time, 100*twang_time/total_time))
  cat(sprintf("    Kbal:           %.2f sec (%.1f%%)\n",
              kbal_time, 100*kbal_time/total_time))
  cat("----------------------------------------\n")

  # ========================================================================
  # ASSEMBLE RESULTS
  # ========================================================================
  res <- tibble(
    runid = i,
    seed = 1000 + i,
    elapsed_time_secs = total_time,
    time_data = data_time,
    time_cf = cf_time,
    time_twang = twang_time,
    time_kbal = kbal_time,
    nc = nc,
    nt = nt,
    f0_sd = f0_sd,
    num_bins = num_bins,
    true_ATT = true_ATT,
    causal_forest = att_cf,
    twang = att_twang,
    kbal = att_kbal
  )

  # Save iteration results
  FNAME_ITER <- here::here(
    "data", "outputs", "sim_toy_results", "toy_new_methods_seq",
    sprintf("run_%03d.csv", i)
  )
  write_csv(res, FNAME_ITER)

  all_results[[i]] <- res
}

# ========================================================================
# COMBINE ALL ITERATION FILES
# ========================================================================
cat("\n========================================\n")
cat("Combining all iteration files...\n")
cat("========================================\n")

combined_results <- bind_rows(all_results) %>%
  write_csv("data/outputs/sim_toy_results/toy_new_methods_seq.csv")

cat("\n=== New Methods Simulation Complete ===\n")
cat(sprintf("Total iterations: %d\n", n_iter))
cat(sprintf("Total time: %.2f minutes\n", sum(combined_results$elapsed_time_secs)/60))
cat("\nTiming summary per iteration:\n")
timing_summary <- combined_results %>%
  summarize(
    mean_total = mean(elapsed_time_secs),
    mean_cf = mean(time_cf),
    mean_twang = mean(time_twang),
    mean_kbal = mean(time_kbal),
    sd_total = sd(elapsed_time_secs)
  )
print(timing_summary)

cat("\nMethod performance summary:\n")
performance_summary <- combined_results %>%
  summarize(
    cf_rmse = sqrt(mean((causal_forest - true_ATT)^2, na.rm = TRUE)),
    cf_bias = mean(causal_forest - true_ATT, na.rm = TRUE),
    cf_na = sum(is.na(causal_forest)),
    twang_rmse = sqrt(mean((twang - true_ATT)^2, na.rm = TRUE)),
    twang_bias = mean(twang - true_ATT, na.rm = TRUE),
    twang_na = sum(is.na(twang)),
    kbal_rmse = sqrt(mean((kbal - true_ATT)^2, na.rm = TRUE)),
    kbal_bias = mean(kbal - true_ATT, na.rm = TRUE),
    kbal_na = sum(is.na(kbal))
  )
print(performance_summary)
