# Simulation script: scripts/sims-bias_mse/sim_toy_new_methods.R
# Run the three new methods on the same seeds as the main simulation

require(tidyverse)
require(mvtnorm)
require(grf)
require(twang)
require(kbal)
require(tictoc)
library(CSM)

# USE_PARALLEL = F
USE_PARALLEL = T
if (USE_PARALLEL) {
  library(foreach)
  library(doParallel)
  cores <- detectCores()
  num_cores <- min(cores - 1, 20)

  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
}

# Create output directory
dir.create("data/outputs/sim_toy_results/toy_new_methods",
           showWarnings = FALSE, recursive = TRUE)

source(here::here("scripts/wrappers.R"))
source(here::here("R/utils.R"))

# Run simulations for new methods only
res <- foreach(
  i = 1:60,
  .packages = c("tidyverse", "mvtnorm", "grf", "twang", "kbal",
                "tictoc", "here", "CSM"),
  .combine = rbind
) %dopar% {

  # Use same seed as main simulation
  set.seed(1000 + i)

  start_time <- Sys.time()
  tic()

  # ========================================================================
  # GENERATE DATA (same as main simulation)
  # ========================================================================
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

  my_covs <- c("X1", "X2")
  zform2 <- as.formula("Z ~ X1+X2")

  # ========================================================================
  # COMPUTE TRUE ATT
  # ========================================================================
  true_ATT <- df %>%
    filter(Z) %>%
    summarize(att = mean(Y1 - Y0)) %>%
    pull(att)

  # ========================================================================
  # COMPUTE ESTIMATES - NEW METHODS
  # ========================================================================
  cat(sprintf("Run %d: Computing Causal Forest...\n", i))
  att_cf <- tryCatch(
    get_att_causal_forest(df, covs = my_covs),
    error = function(e) {
      warning(sprintf("Causal Forest failed in run %d: %s", i, e$message))
      return(NA_real_)
    }
  )

  cat(sprintf("Run %d: Computing TWANG...\n", i))
  att_twang <- tryCatch(
    get_att_twang(df, form = zform2),
    error = function(e) {
      warning(sprintf("TWANG failed in run %d: %s", i, e$message))
      return(NA_real_)
    }
  )

  cat(sprintf("Run %d: Computing Kbal...\n", i))
  att_kbal <- tryCatch(
    get_att_kbal(df, covs = my_covs),
    error = function(e) {
      warning(sprintf("Kbal failed in run %d: %s", i, e$message))
      return(NA_real_)
    }
  )

  # ========================================================================
  # CALCULATE TIMING
  # ========================================================================
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # ========================================================================
  # ASSEMBLE RESULTS
  # ========================================================================
  cat(sprintf("Run %d: Assembling results...\n", i))

  res <- tibble(
    runid = i,
    seed = 1000 + i,
    elapsed_time_secs = elapsed_time,
    nc = nc,
    nt = nt,
    f0_sd = f0_sd,
    num_bins = num_bins,
    true_ATT = true_ATT,
    causal_forest = att_cf,
    twang = att_twang,
    kbal = att_kbal
  )

  toc()

  # Save iteration results
  FNAME_ITER <- here::here(
    "data", "outputs", "sim_toy_results", "toy_new_methods",
    sprintf("run_%03d.csv", i)
  )
  write_csv(res, FNAME_ITER)

  cat(sprintf("Run %d: Complete!\n", i))

  res
}

if (USE_PARALLEL) {
  stopCluster(cl)
}

# ========================================================================
# COMBINE ALL ITERATION FILES
# ========================================================================
cat("Combining all iteration files...\n")

combined_results <- list.files(
  "data/outputs/sim_toy_results/toy_new_methods/",
  pattern = "run_.*\\.csv",
  full.names = TRUE
) %>%
  map_df(read_csv) %>%
  write_csv("data/outputs/sim_toy_results/toy_new_methods.csv")

cat("\n=== New Methods Simulation Complete ===\n")
cat("Results summary:\n")
print(summary(combined_results))

# Print timing info
cat("\nTiming summary:\n")
cat("Mean time per iteration:", mean(combined_results$elapsed_time_secs), "seconds\n")
cat("Total time:", sum(combined_results$elapsed_time_secs), "seconds\n")
