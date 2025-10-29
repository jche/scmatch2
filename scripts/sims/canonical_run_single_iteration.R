#!/usr/bin/env Rscript
# scripts/sims/run_single_iteration.R
# Run a single iteration of the comprehensive simulation

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript run_single_iteration.R <sim_type> <iteration_id>

  sim_type: one of 'acic', 'hainmueller', 'kang'
  iteration_id: integer from 1 to n_iter")
}

sim_type <- args[1]
iteration_id <- as.integer(args[2])

cat(sprintf("\n=== Running %s simulation, iteration %d ===\n", sim_type, iteration_id))

# Load required packages
suppressPackageStartupMessages({
  require(tidyverse)
  require(mvtnorm)
  require(optweight)
  require(dbarts)
  require(tmle)
  require(AIPW)
  require(grf)
  require(twang)
  require(kbal)
  require(tictoc)
  require(CSM)
})

# Load helper functions
source(here::here("R/wrappers.R"))
source(here::here("R/utils.R"))
source(here::here("R/sim_data.R"))

# Set superlearner libraries
SL.library1 <- c("SL.mean", "SL.lm", "SL.glm")
SL.library2 <- c("SL.glm", "SL.glmnet", "SL.randomForest", "SL.xgboost")
SL.library3Q <- c("SL.glm", "tmle.SL.dbarts2", "SL.glmnet")
SL.library3g <- c("SL.glm", "tmle.SL.dbarts.k.5", "SL.gam")

# Set up output directory
output_dir <- file.path("data", "outputs", "sim_slurm", sim_type)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ========================================================================
# SET SEED
# ========================================================================
seed_start <- 1000

if (sim_type == "acic") {
  # ACIC uses random seeds from a specific range
  set.seed(seed_start + iteration_id)
  current_seed <- sample(1:100000, 1)
  set.seed(current_seed)
} else {
  current_seed <- seed_start + iteration_id
  set.seed(current_seed)
}

iter_start <- Sys.time()

# ========================================================================
# GENERATE DATA
# ========================================================================
cat("Generating data...\n")
data_start <- Sys.time()

if (sim_type == "acic") {
  n <- 1000
  p <- 10
  model.trt <- "step"
  model.rsp <- "step"

  df <- gen_df_acic(
    model.trt = model.trt,
    root.trt = 0.35,
    overlap.trt = "full",
    model.rsp = model.rsp,
    alignment = 0.75,
    te.hetero = "high",
    random.seed = current_seed,
    n = n,
    p = p
  )

  covs <- df %>%
    select(starts_with("X")) %>%
    colnames()

  zform1 <- as.formula(paste0("Z ~ ", paste0(covs, collapse = "+")))
  zform2 <- as.formula(paste0("Z ~ (", paste0(covs, collapse = "+"), ")^2"))
  form1 <- as.formula(paste0("Y ~ ", paste0(covs, collapse = "+")))
  form2 <- as.formula(paste0("Y ~ (", paste0(covs, collapse = "+"), ")^2"))
  nbins <- 5

} else if (sim_type == "hainmueller") {
  nc <- 250
  nt <- 50

  df <- gen_df_hain(
    nc = nc,
    nt = nt,
    sigma_e = "n100",
    outcome = "nl2",
    sigma_y = 1,
    ATE = 0
  )

  covs <- c("X1", "X2", "X3", "X4", "X5", "X6")
  zform1 <- as.formula("Z ~ X1+X2+X3+X4+X5+X6")
  zform2 <- as.formula("Z ~ (X1+X2+X3+X4+X5+X6)^2")
  form1 <- as.formula("Y ~ X1+X2+X3+X4+X5+X6")
  form2 <- as.formula("Y ~ (X1+X2+X3+X4+X5+X6)^2")
  nbins <- 5

} else if (sim_type == "kang") {
  n <- 1000

  df <- gen_df_kang(n = n)

  covs <- df %>%
    select(starts_with("X")) %>%
    colnames()

  zform1 <- as.formula(paste0("Z ~ ", paste0(covs, collapse = "+")))
  zform2 <- as.formula(paste0("Z ~ (", paste0(covs, collapse = "+"), ")^2"))
  form1 <- as.formula(paste0("Y ~ ", paste0(covs, collapse = "+")))
  form2 <- as.formula(paste0("Y ~ (", paste0(covs, collapse = "+"), ")^2"))
  nbins <- 5

} else {
  stop("sim_type must be one of: 'acic', 'hainmueller', 'kang'")
}

# Common processing
dist_scaling <- df %>%
  summarize(across(
    starts_with("X"),
    function(x) {
      if (is.numeric(x)) nbins / (max(x) - min(x))
      else 1000
    }
  ))

preds_csm <- get_cal_matches(
  df = df,
  metric = "maximum",
  scaling = dist_scaling,
  rad_method = "fixed",
  est_method = "average",
  k = 25
)

preds_cem <- get_cem_matches(
  df = df,
  num_bins = nbins,
  est_method = "average",
  return = "all"
)

ninf <- length(attr(preds_csm, "unmatched_units"))
ninf_cem <- sum(df$Z) - length(attr(preds_cem, "feasible_units"))

df <- df %>%
  filter(!id %in% attr(preds_csm, "unmatched_units"))

data_time <- as.numeric(difftime(Sys.time(), data_start, units = "secs"))

# True ATT
if (sim_type == "acic") {
  true_ATT <- df %>%
    filter(Z) %>%
    summarize(att = mean(Y1 - Y0)) %>%
    pull(att)
} else {
  true_ATT <- 0
}

cat(sprintf("Data generated: n=%d, true_ATT=%.3f, time=%.2fs\n",
            nrow(df), true_ATT, data_time))

# ========================================================================
# COMPUTE ALL ESTIMATES
# ========================================================================

safe_compute <- function(method_name, expr) {
  t_start <- Sys.time()
  result <- tryCatch(
    expr,
    error = function(e) {
      cat(sprintf("  %s: ERROR - %s\n", method_name, e$message))
      return(NA_real_)
    }
  )
  t_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  if (!is.na(result)) {
    cat(sprintf("  %s: %.3f (%.2fs)\n", method_name, result, t_elapsed))
  }
  return(result)
}

cat("\nComputing estimates...\n")

att_diff <- safe_compute("diff", get_att_diff(df))
att_bal1 <- safe_compute("bal1", get_att_bal(df, zform1, rep(0.01, length(covs))))
att_bal2 <- safe_compute("bal2", get_att_bal(df, zform2, rep(0.1, length(covs) + choose(length(covs), 2))))

att_or_lm <- safe_compute("or_lm", get_att_or_lm(df, form = form2))
att_or_bart <- safe_compute("or_bart", get_att_or_bart(df, covs = covs))

att_ps_lm <- safe_compute("ps_lm", get_att_ps_lm(df, zform2))
att_ps_bart <- safe_compute("ps_bart", get_att_ps_bart(df, covs = covs))

att_csm_scm <- safe_compute("csm_scm", get_att_csm(df, scaling = dist_scaling, est_method = "scm", rad_method = "fixed"))
att_csm_avg <- safe_compute("csm_avg", get_att_csm(df, scaling = dist_scaling, est_method = "average", rad_method = "fixed"))
att_cem_scm <- safe_compute("cem_scm", get_att_cem(df, num_bins = nbins, est_method = "scm", estimand = "CEM-ATT"))
att_cem_avg <- safe_compute("cem_avg", get_att_cem(df, num_bins = nbins, est_method = "average", estimand = "CEM-ATT"))
att_onenn <- safe_compute("onenn", get_att_1nn(df, scaling = dist_scaling))

att_tmle1 <- safe_compute("tmle1", get_att_tmle(df, covs = covs, Q.SL.library = SL.library1, g.SL.library = SL.library1))
att_aipw1 <- safe_compute("aipw1", get_att_aipw(df, covs = covs, Q.SL.library = SL.library1, g.SL.library = SL.library1))
att_tmle2 <- safe_compute("tmle2", get_att_tmle(df, covs = covs, Q.SL.library = SL.library3Q, g.SL.library = SL.library3g))
att_aipw2 <- safe_compute("aipw2", get_att_aipw(df, covs = covs, Q.SL.library = SL.library2, g.SL.library = SL.library2))

# New methods
att_cf <- safe_compute("causal_forest", get_att_causal_forest(df, covs = covs))
att_twang <- safe_compute("twang", get_att_twang(df, form = zform1))
att_kbal <- safe_compute("kbal", get_att_kbal(df, covs = covs))

# ========================================================================
# SAVE RESULTS
# ========================================================================
total_time <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))

res <- tibble(
  runid = iteration_id,
  seed = current_seed,
  sim_type = sim_type,
  elapsed_time_secs = total_time,
  ninf = ninf,
  ninf_cem = ninf_cem,
  true_ATT = true_ATT,
  diff = att_diff,
  bal1 = att_bal1,
  bal2 = att_bal2,
  or_lm = att_or_lm,
  or_bart = att_or_bart,
  ps_lm = att_ps_lm,
  ps_bart = att_ps_bart,
  csm_scm = att_csm_scm,
  csm_avg = att_csm_avg,
  cem_scm = att_cem_scm,
  cem_avg = att_cem_avg,
  onenn = att_onenn,
  tmle1 = att_tmle1,
  aipw1 = att_aipw1,
  tmle2 = att_tmle2,
  aipw2 = att_aipw2,
  causal_forest = att_cf,
  twang = att_twang,
  kbal = att_kbal
)

output_file <- file.path(output_dir, sprintf("iter_%04d.csv", iteration_id))
write_csv(res, output_file)

cat(sprintf("\nIteration %d complete! Total time: %.2f seconds\n", iteration_id, total_time))
cat(sprintf("Results saved to: %s\n", output_file))
