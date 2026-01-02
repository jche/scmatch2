#!/usr/bin/env Rscript

# scripts/sims-variance/4_debug_one_iter_V_vs_VE.R
#
# Usage:
#   Rscript scripts/sims-variance/4_debug_one_iter_V_vs_VE.R 1 very_low homoskedastic
#
# What it does:
#   - Runs ONE toy iteration (same generator + matcher)
#   - Runs TWO matching procedures side-by-side:
#       (A) adaptive radius + SCM weights (your current default)
#       (B) k-NN (k=5) + average weights
#   - Computes V_E (pooled) via get_measurement_error_variance()
#   - Manually reconstructs the pieces used inside get_total_variance:
#       squared_deviations, control_weight_correction, V_correction, V_total
#   - Prints scale checks + SE comparisons for both methods in a single table
#   - Prints reuse diagnostics + subclass control counts for both methods
#
# Notes:
#   - Your "true variance" check: if f0_sd is an SD, then true sigma^2 is f0_sd^2.
#     This script prints that implied value for clarity.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(here)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

# ---- Parse args ----
# args <- commandArgs(trailingOnly = TRUE)
args <- NULL
iter <- as.integer(args[[1]] %||% 1)
overlap_label <- args[[2]] %||% "very_low"
error_label <- args[[3]] %||% "homoskedastic"

# ---- Load package code + sim utils ----
devtools::load_all()
source(here::here("scripts/sims-variance/0_sim_inference_utils.R"))

# --- seed mapping (copied from sim_master) ---
generate_seed <- function(i, overlap_label, error_label, N = 600) {
  overlap_levels <- c("very_low", "low", "mid", "high", "very_high")
  error_levels <- c("homoskedastic", "covariate_dep", "treatment_dep")
  N_levels <- c(600)

  overlap_index <- match(overlap_label, overlap_levels); if (is.na(overlap_index)) overlap_index <- 1
  error_index <- match(error_label, error_levels); if (is.na(error_index)) error_index <- 1
  N_index <- match(N, N_levels); if (is.na(N_index)) N_index <- 1

  i + 1000 * (overlap_index - 1) + 10000 * (error_index - 1) + 100000 * (N_index - 1)
}

seed <- generate_seed(iter, overlap_label, error_label, N = 600)

prop_nc_unif_values <- c(
  very_low   = 2/3,
  low        = 1/2,
  mid        = 1/3,
  high       = 1/5,
  very_high  = 1/10
)
prop_nc_unif <- prop_nc_unif_values[[overlap_label]]
stopifnot(!is.null(prop_nc_unif))

# ---- Pretty header ----
cat("\n============================================================\n")
cat("DEBUG ONE ITERATION (COMPARE MATCHING METHODS)\n")
cat("  iter         =", iter, "\n")
cat("  overlap      =", overlap_label, "(prop_nc_unif=", prop_nc_unif, ")\n")
cat("  error_label  =", error_label, "\n")
cat("  seed         =", seed, "\n")
cat("============================================================\n\n")

# ---- Generate df (same as toy_match_infer) ----
f0_sd <- 0.5

df <- make_csm_toy_df(
  nc = 500, nt = 100,
  f0_sd = f0_sd,
  prop_nc_unif = prop_nc_unif,
  ctr_dist = 0.5,
  seed = seed
)

cat("---- Generator variance sanity check ----\n")
cat("f0_sd =", f0_sd, " => implied sigma^2 =", sprintf("%.6f", f0_sd^2), "\n\n")

scaling <- compute_toy_scaling(df, nbins = 6)

# ------------------------------------------------------------------
# Helper: compute all debugging stats from matches_df
# ------------------------------------------------------------------
compute_debug_stats <- function(matches_df, outcome = "Y", treatment = "Z") {

  # core point estimate
  att_hat <- get_att_point_est(matches_df, treatment = treatment, outcome = outcome)

  # pooled V_E as implemented
  ve_res <- get_measurement_error_variance(
    matches_table = matches_df,
    outcome = outcome,
    treatment = treatment,
    var_weight_type = "ess_units"
  )

  V_E <- ve_res$V_E
  sigma2_hat <- ve_res$sigma_hat^2
  N_T <- ve_res$N_T
  ESS_C <- ve_res$ESS_C

  # individual effects (same path as get_total_variance)
  tx_units <- matches_df %>%
    filter(.data[[treatment]] == 1) %>%
    select(id, subclass, .data[[outcome]]) %>%
    rename(Y_t = .data[[outcome]])

  control_outcomes <- matches_df %>%
    filter(.data[[treatment]] == 0) %>%
    group_by(subclass) %>%
    summarize(Y_hat_0 = sum(.data[[outcome]] * weights) / sum(weights), .groups = "drop")

  individual_effects <- tx_units %>%
    left_join(control_outcomes, by = "subclass") %>%
    mutate(indiv_effect = Y_t - Y_hat_0)

  # squared deviations term (as in get_total_variance)
  squared_deviations <- sum((individual_effects$indiv_effect - att_hat)^2) / N_T
  var_indiv <- var(individual_effects$indiv_effect)

  # correction term (as in get_total_variance)
  control_weight_correction <- matches_df %>%
    filter(.data[[treatment]] == 0) %>%
    group_by(id) %>%
    summarize(
      sum_weights = sum(weights),
      sum_squared_weights = sum(weights^2),
      .groups = "drop"
    ) %>%
    summarize(
      correction_term = sum((sum_weights^2 - sum_squared_weights)) / N_T
    ) %>%
    pull(correction_term)

  V_correction <- sigma2_hat * control_weight_correction
  V_total <- squared_deviations + V_correction

  # scale checks
  NT_VE <- N_T * V_E
  ratio_Vcorr_to_NT_VE <- V_correction / NT_VE
  ratio_Vtot_to_NT_VE  <- V_total / NT_VE

  # SE diagnostics
  SE_code   <- sqrt(V_total) / sqrt(N_T)          # matches your get_total_variance return scaling
  SE_VEonly <- sqrt(V_E)                          # your pooled_V_E_CI approach
  SE_alt1   <- sqrt(squared_deviations / N_T + V_E)
  SE_alt2   <- sqrt(var_indiv / N_T + V_E)

  tibble(
    ATT_hat = att_hat,
    N_T = N_T,
    ESS_C = ESS_C,
    sigma2_hat = sigma2_hat,
    V_E = V_E,
    squared_deviations = squared_deviations,
    var_indiv = var_indiv,
    control_weight_correction = control_weight_correction,
    V_correction = V_correction,
    V_total = V_total,
    NT_times_VE = NT_VE,
    ratio_Vcorr_to_NT_VE = ratio_Vcorr_to_NT_VE,
    ratio_Vtot_to_NT_VE  = ratio_Vtot_to_NT_VE,
    SE_code = SE_code,
    SE_VEonly = SE_VEonly,
    SE_alt1 = SE_alt1,
    SE_alt2 = SE_alt2
  )
}

# ------------------------------------------------------------------
# Helper: reuse + small subclass counts diagnostics
# ------------------------------------------------------------------
compute_extra_diags <- function(matches_df, treatment = "Z") {

  reuse_tbl <- matches_df %>%
    filter(.data[[treatment]] == 0) %>%
    group_by(id) %>%
    summarise(sum_w = sum(weights), .groups = "drop")

  subclass_ctrl_n <- matches_df %>%
    filter(.data[[treatment]] == 0) %>%
    count(subclass, name = "n_ctrl") %>%
    arrange(n_ctrl)

  list(
    n_unique_controls_used = nrow(reuse_tbl),
    reuse_summary = summary(reuse_tbl$sum_w),
    reuse_top10 = reuse_tbl %>% arrange(desc(sum_w)) %>% head(10),
    subclass_ctrl_n_head20 = head(subclass_ctrl_n, 20)
  )
}

# ------------------------------------------------------------------
# Run TWO matchers
# ------------------------------------------------------------------

# (A) adaptive + SCM (your current config)
mtch_scm <- get_cal_matches(
  data = df,
  form = Z ~ X1 + X2,
  rad_method = "adaptive",
  scaling = scaling,
  k = 5,
  warn = FALSE,
  est_method = "scm"
)
matches_scm <- full_unit_table(mtch_scm)

# (B) 5-NN + average weights
mtch_avg <- get_cal_matches(
  data = df,
  form = Z ~ X1 + X2,
  rad_method = "knn",
  scaling = scaling,
  k = 5,
  warn = FALSE,
  est_method = "average"
)
matches_avg <- full_unit_table(mtch_avg)

# ------------------------------------------------------------------
# Compute and print comparison table
# ------------------------------------------------------------------
stats_scm <- compute_debug_stats(matches_scm) %>% mutate(method = "adaptive + scm")
stats_avg <- compute_debug_stats(matches_avg) %>% mutate(method = "knn(k=5) + average")

compare_tbl <- bind_rows(stats_scm, stats_avg) %>%
  select(method, everything())

cat("---- Summary comparison (both methods) ----\n")
print(compare_tbl)

# ------------------------------------------------------------------
# Print the same info in your old narrative style (per method)
# ------------------------------------------------------------------
print_method_block <- function(method_name, stats_row) {
  cat("\n============================================================\n")
  cat("METHOD:", method_name, "\n")
  cat("============================================================\n")

  cat("---- Point estimate ----\n")
  cat("ATT_hat =", sprintf("%.6f", stats_row$ATT_hat), "\n\n")

  cat("---- Sample sizes ----\n")
  cat("N_T   =", stats_row$N_T, "\n")
  cat("ESS_C =", sprintf("%.3f", stats_row$ESS_C), "\n\n")

  cat("---- Measurement error piece (pooled) ----\n")
  cat("sigma_hat^2 =", sprintf("%.6f", stats_row$sigma2_hat), "\n")
  cat("V_E         =", sprintf("%.8f", stats_row$V_E), "  (ATT scale)\n\n")

  cat("---- Heterogeneity / empirical piece ----\n")
  cat("squared_deviations (mean((tau_i-att)^2)) =", sprintf("%.8f", stats_row$squared_deviations), "\n")
  cat("var(indiv_effect)                        =", sprintf("%.8f", stats_row$var_indiv), "\n")
  cat("check squared_deviations / var           =", sprintf("%.6f", stats_row$squared_deviations / stats_row$var_indiv), "\n\n")

  cat("---- Weight correction piece ----\n")
  cat("control_weight_correction =", sprintf("%.8f", stats_row$control_weight_correction), "\n")
  cat("V_correction              =", sprintf("%.8f", stats_row$V_correction), "\n\n")

  cat("---- Total V (as coded) ----\n")
  cat("V_total =", sprintf("%.8f", stats_row$V_total), "\n\n")

  cat("---- Unit / double-count checks ----\n")
  cat("N_T * V_E                  =", sprintf("%.8f", stats_row$NT_times_VE), "\n")
  cat("V_correction / (N_T * V_E)  =", sprintf("%.6f", stats_row$ratio_Vcorr_to_NT_VE), "\n")
  cat("V_total / (N_T * V_E)       =", sprintf("%.6f", stats_row$ratio_Vtot_to_NT_VE), "\n\n")

  cat("---- SE comparisons (diagnostic) ----\n")
  cat("SE_code    = sqrt(V_total)/sqrt(N_T)                 =", sprintf("%.6f", stats_row$SE_code), "\n")
  cat("SE_VEonly  = sqrt(V_E)                               =", sprintf("%.6f", stats_row$SE_VEonly), "\n")
  cat("SE_alt1    = sqrt(squared_deviations/N_T + V_E)      =", sprintf("%.6f", stats_row$SE_alt1), "\n")
  cat("SE_alt2    = sqrt(var(indiv_effect)/N_T + V_E)       =", sprintf("%.6f", stats_row$SE_alt2), "\n")
}

print_method_block("adaptive + scm", stats_scm %>% select(-method) %>% as.list())
print_method_block("knn(k=5) + average", stats_avg %>% select(-method) %>% as.list())

# ------------------------------------------------------------------
# Extra diagnostics (reuse + subclass counts)
# ------------------------------------------------------------------
diag_scm <- compute_extra_diags(matches_scm)
diag_avg <- compute_extra_diags(matches_avg)

cat("\n============================================================\n")
cat("EXTRA DIAGNOSTICS: CONTROL REUSE\n")
cat("============================================================\n")

cat("\n---- adaptive + scm ----\n")
cat("controls used (unique ids) =", diag_scm$n_unique_controls_used, "\n")
cat("summary(sum_w):\n")
print(diag_scm$reuse_summary)
cat("\nTop 10 reused controls by sum_w:\n")
print(diag_scm$reuse_top10)

cat("\n---- knn(k=5) + average ----\n")
cat("controls used (unique ids) =", diag_avg$n_unique_controls_used, "\n")
cat("summary(sum_w):\n")
print(diag_avg$reuse_summary)
cat("\nTop 10 reused controls by sum_w:\n")
print(diag_avg$reuse_top10)

cat("\n============================================================\n")
cat("EXTRA DIAGNOSTICS: CONTROL COUNTS PER SUBCLASS (SMALLEST FIRST)\n")
cat("============================================================\n")

cat("\n---- adaptive + scm (head 20) ----\n")
print(diag_scm$subclass_ctrl_n_head20)

cat("\n---- knn(k=5) + average (head 20) ----\n")
print(diag_avg$subclass_ctrl_n_head20)

cat("\nDONE.\n")
