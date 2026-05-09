#!/usr/bin/env Rscript
# scripts/sims-variance/4_make_latex_tables.R
#
# Produces LaTeX for Tables 1, 3-8 in variance_estimation_memo.tex.
#
# Prerequisites: run 2_collect_results.R for each dataset first.
#   Rscript scripts/sims-variance/2_collect_results.R sims-variance
#   Rscript scripts/sims-variance/2_collect_results.R "sims-variance-fully-adaptive(1nn)"
#   Rscript scripts/sims-variance/2_collect_results.R "sims-variance-fully-adaptive(2nn)"
#
# Usage:
#   Rscript scripts/sims-variance/4_make_latex_tables.R
#
# Output: prints LaTeX blocks to stdout, one per table.
#
# Table mapping:
#   Table 1  tab:toy_var_pooled_vs_het      k=5, PATT (vs ATT=3)
#   Table 3  tab:toy_var_pooled_vs_het-V_E  k=5, SATT (V_E-only CI)
#   Table 4  tab:toy_var_pooled_vs_het_1nn  k=1, PATT
#   Table 5  tab:toy_var_pooled_vs_het_VE_1nn  k=1, SATT
#   Table 6  tab:toy_var_pooled_vs_het_2nn  k=2, PATT
#   Table 7  tab:toy_var_pooled_vs_het_VE_2nn  k=2, SATT
#   Table 8  tab:compare_k125              k=1/2/5 SATT comparison

suppressPackageStartupMessages({
  library(dplyr)
  library(here)
})

ATT_POP_TRUTH <- 3

overlap_levels <- c("very_low", "low", "mid", "high", "very_high")
overlap_labels <- c(
  very_low  = "Very Low",
  low       = "Low",
  mid       = "Medium",
  high      = "High",
  very_high = "Very High"
)

# -------------------------------------------------------------------
# Load datasets
# -------------------------------------------------------------------
load_combined <- function(output_name) {
  rds <- here::here("data/outputs", output_name, "combined_results.rds")
  if (!file.exists(rds))
    stop("Not found: ", rds, "\nRun: Rscript scripts/sims-variance/2_collect_results.R \"", output_name, "\"")
  df <- readRDS(rds)
  df$deg_overlap <- factor(df$deg_overlap, levels = overlap_levels)
  df
}

message("Loading sims-variance (k=5) ...")
df5 <- load_combined("sims-variance")

message("Loading sims-variance-fully-adaptive(1nn) (k=1) ...")
df1 <- load_combined("sims-variance-fully-adaptive(1nn)")

message("Loading sims-variance-fully-adaptive(2nn) (k=2) ...")
df2 <- load_combined("sims-variance-fully-adaptive(2nn)")

# -------------------------------------------------------------------
# Summarise helpers
# -------------------------------------------------------------------

# PATT table: methods pooled (homo) and pooled_het (het), CI vs ATT=3
summarise_patt <- function(df) {
  df %>%
    filter(inference_method %in% c("pooled", "pooled_het"),
           !is.na(att_est), !is.na(SE)) %>%
    mutate(
      covers_pop = CI_lower <= ATT_POP_TRUTH & ATT_POP_TRUTH <= CI_upper,
      error_pop  = att_est - ATT_POP_TRUTH
    ) %>%
    group_by(deg_overlap, inference_method) %>%
    summarise(
      N_C_tilde  = mean(ESS_C,      na.rm = TRUE),
      true_SE    = sd(error_pop,    na.rm = TRUE),
      avg_SE_hat = mean(SE,         na.rm = TRUE),
      coverage   = mean(covers_pop, na.rm = TRUE),
      n_runs     = sum(!is.na(att_est)),
      .groups    = "drop"
    )
}

# SATT table: methods pooled_V_E_CI (homo) and pooled_het_V_E_CI (het), CI vs SATT
summarise_satt <- function(df) {
  df %>%
    filter(inference_method %in% c("pooled_V_E_CI", "pooled_het_V_E_CI"),
           !is.na(att_est), !is.na(SE), !is.na(SATT)) %>%
    mutate(
      covers_satt = CI_lower <= SATT & SATT <= CI_upper,
      error_satt  = att_est - SATT
    ) %>%
    group_by(deg_overlap, inference_method) %>%
    summarise(
      N_C_tilde  = mean(ESS_C,       na.rm = TRUE),
      true_SE    = sd(error_satt,    na.rm = TRUE),
      avg_SE_hat = mean(SE,          na.rm = TRUE),
      coverage   = mean(covers_satt, na.rm = TRUE),
      n_runs     = sum(!is.na(att_est)),
      .groups    = "drop"
    )
}

# -------------------------------------------------------------------
# LaTeX printer for standard 6-column tables
# (Overlap | N_C_tilde | True SE | Homo Est.SE | Homo Cov | Het Est.SE | Het Cov)
# -------------------------------------------------------------------
print_standard_table <- function(summ, label, caption, homo_method, het_method) {
  # Sort by descending N_C_tilde so rows are always labeled
  # Very High (largest control pool) ... Very Low (smallest), regardless of
  # whether deg_overlap labels are inverted in the raw data (k=5 dataset has
  # inverted labels: very_high -> N_C~17, very_low -> N_C~53).
  display_labels <- c("Very High", "High", "Medium", "Low", "Very Low")

  homo_rows <- summ %>% filter(inference_method == homo_method) %>% arrange(desc(N_C_tilde))
  het_rows  <- summ %>% filter(inference_method == het_method)  %>% arrange(desc(N_C_tilde))

  fmt3 <- function(x) if (is.na(x)) "---" else sprintf("%.3f", x)
  fmt1 <- function(x) if (is.na(x)) "---" else sprintf("%.1f", x)

  cat(sprintf("%% Label: %s\n", label))
  cat(sprintf("%% Run: Rscript scripts/sims-variance/4_make_latex_tables.R\n\n"))
  cat("\\begin{table}[hbt!]\n")
  cat("\\centering\n")
  cat("\\begin{tabular}{lcccccc}\n")
  cat("  \\hline\n")
  cat("  & & & \\multicolumn{2}{c}{Homoskedastic} & \\multicolumn{2}{c}{Heteroskedastic} \\\\\n")
  cat("  \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}\n")
  cat("  Overlap & $\\tilde{N}_C$ & True SE & Est. SE & Coverage & Est. SE & Coverage \\\\\n")
  cat("  \\hline\n")

  for (i in seq_len(nrow(homo_rows))) {
    lab      <- display_labels[i]
    nc       <- homo_rows$N_C_tilde[i]
    true_se  <- homo_rows$true_SE[i]
    homo_se  <- homo_rows$avg_SE_hat[i]
    homo_cov <- homo_rows$coverage[i]
    het_se   <- het_rows$avg_SE_hat[i]
    het_cov  <- het_rows$coverage[i]

    cat(sprintf("  %s & %s & %s & %s & %s & %s & %s \\\\\n",
                lab, fmt1(nc), fmt3(true_se),
                fmt3(homo_se), fmt3(homo_cov),
                fmt3(het_se),  fmt3(het_cov)))
  }

  cat("   \\hline\n")
  cat("\\end{tabular}\n")
  cat(sprintf("\\caption{%s}\n", caption))
  cat(sprintf("\\label{%s}\n", label))
  cat("\\end{table}\n\n")
}

# -------------------------------------------------------------------
# Table 1: PATT, k=5
# -------------------------------------------------------------------
n5 <- length(unique(df5$runID))
summ5_patt <- summarise_patt(df5)

sep <- paste(rep("-", 60), collapse = "")
cat("\n% ", sep, "\n")
cat("% TABLE 1 (tab:toy_var_pooled_vs_het): PATT estimators, k=5 adaptive\n")
cat("% ", sep, "\n\n")

print_standard_table(
  summ5_patt,
  label      = "tab:toy_var_pooled_vs_het",
  caption    = sprintf(paste0(
    "Performance of variance estimators under varying degrees of overlap. ",
    "$\\tilde{N}_C$ denotes the average effective sample size of the control group. ",
    "True SE is the Monte Carlo standard deviation of $\\hat{\\tau}$ across %d simulations, ",
    "and Est. SE is the average of the estimated standard errors. ",
    "``Homoskedastic'' pools residual variance across matched sets; ",
    "``Heteroskedastic'' applies the heteroskedasticity-consistent correction."), n5),
  homo_method = "pooled",
  het_method  = "pooled_het"
)

# -------------------------------------------------------------------
# Table 3: SATT (V_E-only), k=5
# -------------------------------------------------------------------
summ5_satt <- summarise_satt(df5)

cat("\n% ", sep, "\n")
cat("% TABLE 3 (tab:toy_var_pooled_vs_het-V_E): SATT V_E estimators, k=5 adaptive\n")
cat("% ", sep, "\n\n")

print_standard_table(
  summ5_satt,
  label      = "tab:toy_var_pooled_vs_het-V_E",
  caption    = sprintf(paste0(
    "Performance of variance estimators under varying degrees of overlap. ",
    "$\\tilde{N}_C$ denotes the average effective sample size of the control group. ",
    "True SE is the Monte Carlo standard deviation of $\\hat{\\tau}$ across %d simulations, ",
    "and Est. SE is the average of the estimated standard errors. ",
    "``Homoskedastic'' pools residual variance across matched sets; ",
    "``Heteroskedastic'' applies the heteroskedasticity-consistent correction."), n5),
  homo_method = "pooled_V_E_CI",
  het_method  = "pooled_het_V_E_CI"
)

# -------------------------------------------------------------------
# Table 4: PATT, k=1
# -------------------------------------------------------------------
n1 <- length(unique(df1$runID))
summ1_patt <- summarise_patt(df1)

cat("\n% ", sep, "\n")
cat("% TABLE 4 (tab:toy_var_pooled_vs_het_1nn): PATT estimators, k=1 adaptive\n")
cat("% ", sep, "\n\n")

print_standard_table(
  summ1_patt,
  label      = "tab:toy_var_pooled_vs_het_1nn",
  caption    = sprintf(paste0(
    "Performance of population-variance estimators under the fully adaptive ",
    "one-nearest-neighbor design. ",
    "$\\tilde{N}_C$ denotes the average effective sample size of the control group. ",
    "True SE is the Monte Carlo standard deviation of $\\hat{\\tau}$ across %d simulations, ",
    "and Est. SE is the average estimated standard error. ",
    "``Homoskedastic'' pools residual variance across matched sets; ",
    "``Heteroskedastic'' applies the heteroskedasticity-consistent correction."), n1),
  homo_method = "pooled",
  het_method  = "pooled_het"
)

# -------------------------------------------------------------------
# Table 5: SATT (V_E-only), k=1
# -------------------------------------------------------------------
summ1_satt <- summarise_satt(df1)

cat("\n% ", sep, "\n")
cat("% TABLE 5 (tab:toy_var_pooled_vs_het_VE_1nn): SATT V_E estimators, k=1 adaptive\n")
cat("% ", sep, "\n\n")

print_standard_table(
  summ1_satt,
  label      = "tab:toy_var_pooled_vs_het_VE_1nn",
  caption    = sprintf(paste0(
    "Performance of SATT-oriented $V_E$ estimators under the fully adaptive ",
    "one-nearest-neighbor design. ",
    "$\\tilde{N}_C$ denotes the average effective sample size of the control group. ",
    "True SE is the Monte Carlo standard deviation of $\\hat{\\tau}$ across %d simulations, ",
    "and Est. SE is the average estimated standard error. ",
    "``Homoskedastic'' pools residual variance across matched sets; ",
    "``Heteroskedastic'' applies the heteroskedasticity-consistent correction."), n1),
  homo_method = "pooled_V_E_CI",
  het_method  = "pooled_het_V_E_CI"
)

# -------------------------------------------------------------------
# Table 6: PATT, k=2
# -------------------------------------------------------------------
n2 <- length(unique(df2$runID))
summ2_patt <- summarise_patt(df2)

cat("\n% ", sep, "\n")
cat("% TABLE 6 (tab:toy_var_pooled_vs_het_2nn): PATT estimators, k=2 adaptive\n")
cat("% ", sep, "\n\n")

print_standard_table(
  summ2_patt,
  label      = "tab:toy_var_pooled_vs_het_2nn",
  caption    = sprintf(paste0(
    "Performance of population-variance (PATT) estimators under the fully adaptive ",
    "two-nearest-neighbor design. ",
    "$\\tilde{N}_C$ is the average effective control sample size. ",
    "True SE is the Monte Carlo SD of $\\hat{\\tau}$ across $\\approx %d$ simulations. ",
    "``Homoskedastic'' pools residual variance; ",
    "``Heteroskedastic'' applies the het-consistent correction."), n2),
  homo_method = "pooled",
  het_method  = "pooled_het"
)

# -------------------------------------------------------------------
# Table 7: SATT (V_E-only), k=2
# -------------------------------------------------------------------
summ2_satt <- summarise_satt(df2)

cat("\n% ", sep, "\n")
cat("% TABLE 7 (tab:toy_var_pooled_vs_het_VE_2nn): SATT V_E estimators, k=2 adaptive\n")
cat("% ", sep, "\n\n")

print_standard_table(
  summ2_satt,
  label      = "tab:toy_var_pooled_vs_het_VE_2nn",
  caption    = sprintf(paste0(
    "Performance of SATT-oriented $\\hat{V}_E$ estimators under the fully adaptive ",
    "two-nearest-neighbor design. ",
    "$\\tilde{N}_C$ is the average effective control sample size. ",
    "True SE is the Monte Carlo SD of $\\hat{\\tau}$ across $\\approx %d$ simulations. ",
    "``Homoskedastic'' pools residual variance; ",
    "``Heteroskedastic'' applies the het-consistent correction."), n2),
  homo_method = "pooled_V_E_CI",
  het_method  = "pooled_het_V_E_CI"
)

# -------------------------------------------------------------------
# Table 8: tab:compare_k125
# Compare het vs pooled SATT V_E coverage across k=1, 2, 5.
# Rows are aligned by similar N_C_tilde (ascending order).
# k=1/2 ordering: very_low < low < mid < high < very_high
# k=5   ordering: very_high < high < mid < low < very_low  (reversed)
# -------------------------------------------------------------------

get_satt_cov <- function(df) {
  df %>%
    filter(inference_method %in% c("pooled_V_E_CI", "pooled_het_V_E_CI"),
           !is.na(att_est), !is.na(SATT)) %>%
    mutate(covers_satt = CI_lower <= SATT & SATT <= CI_upper) %>%
    group_by(deg_overlap, inference_method) %>%
    summarise(
      N_C_tilde = mean(ESS_C,       na.rm = TRUE),
      coverage  = mean(covers_satt, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = inference_method,
                       values_from = coverage,
                       names_prefix = "cov_") %>%
    arrange(N_C_tilde)
}

suppressPackageStartupMessages(library(tidyr))

tab8_k1 <- get_satt_cov(df1)
tab8_k2 <- get_satt_cov(df2)
tab8_k5 <- get_satt_cov(df5)

# All three should have 5 rows (one per overlap level), ordered by N_C_tilde ascending.
# Row i of k=1 pairs with row i of k=2 and row i of k=5.
stopifnot(nrow(tab8_k1) == 5, nrow(tab8_k2) == 5, nrow(tab8_k5) == 5)

fmt3 <- function(x) if (is.na(x)) "---" else sprintf("%.3f", x)
fmtnc <- function(nc1, nc2, nc5) {
  vals <- sort(c(nc1, nc2, nc5))
  lo <- floor(min(vals))
  hi <- ceiling(max(vals))
  if (lo == hi) sprintf("$\\approx %d$", lo) else sprintf("$\\approx %d$--%d", lo, hi)
}

cat("\n% ", sep, "\n")
cat("% TABLE 8 (tab:compare_k125): SATT V_E coverage comparison k=1,2,5\n")
cat("% ", sep, "\n\n")

cat("% Label: tab:compare_k125\n")
cat("% Run: Rscript scripts/sims-variance/4_make_latex_tables.R\n\n")

cat("\\begin{table}[hbt!]\n")
cat("\\centering\n")
cat("\\begin{tabular}{ccccccc}\n")
cat("  \\hline\n")
cat("  & \\multicolumn{2}{c}{$k=1$} & \\multicolumn{2}{c}{$k=2$} & \\multicolumn{2}{c}{$k=5$} \\\\\n")
cat("  \\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\n")
cat("  $\\tilde{N}_C$ & Het Cov & Pooled Cov & Het Cov & Pooled Cov & Het Cov & Pooled Cov \\\\\n")
cat("  \\hline\n")

for (i in 1:5) {
  nc_lbl <- fmtnc(tab8_k1$N_C_tilde[i], tab8_k2$N_C_tilde[i], tab8_k5$N_C_tilde[i])
  het1   <- fmt3(tab8_k1$cov_pooled_het_V_E_CI[i])
  pool1  <- fmt3(tab8_k1$cov_pooled_V_E_CI[i])
  het2   <- fmt3(tab8_k2$cov_pooled_het_V_E_CI[i])
  pool2  <- fmt3(tab8_k2$cov_pooled_V_E_CI[i])
  het5   <- fmt3(tab8_k5$cov_pooled_het_V_E_CI[i])
  pool5  <- fmt3(tab8_k5$cov_pooled_V_E_CI[i])

  cat(sprintf("  %s & %s & %s & %s & %s & %s & %s \\\\\n",
              nc_lbl, het1, pool1, het2, pool2, het5, pool5))
}

cat("  \\hline\n")
cat("\\end{tabular}\n")
caption8 <- paste0(
  "SATT-oriented $\\hat{V}_E$ coverage for the heteroskedastic and pooled estimators ",
  "across comparable effective sample sizes $\\tilde{N}_C$, for $k=1$, $k=2$, and $k=5$ ",
  "adaptive-caliper SCM matching. Each row pairs scenarios from the three designs at similar ",
  "$\\tilde{N}_C$. Nominal coverage is 0.95."
)
cat(sprintf("\\caption{%s}\n", caption8))
cat("\\label{tab:compare_k125}\n")
cat("\\end{table}\n\n")

message("Done. All 7 tables printed.")
