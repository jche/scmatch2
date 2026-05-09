# scripts/sims-variance-het-sigma/3_summarize_results.R
#!/usr/bin/env Rscript
#
# Usage: Rscript 3_summarize_results.R [output_name]

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(here)
})

args        <- commandArgs(trailingOnly = TRUE)
output_name <- if (length(args) >= 1) args[[1]] else "sims-variance-het-sigma"

in_csv  <- here::here(file.path("data/outputs", output_name, "combined_results.csv"))
out_dir <- here::here(file.path("tables", output_name))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

methods <- c("homo", "het", "alt_common", "alt_tt")

overlap_levels <- c("very_low", "low", "mid", "high", "very_high")
overlap_labels <- c(
  very_low  = "Very Low",
  low       = "Low",
  mid       = "Medium",
  high      = "High",
  very_high = "Very High"
)

ATT_POP_TRUTH <- 3

stopifnot(file.exists(in_csv))
res <- readr::read_csv(in_csv, show_col_types = FALSE) %>% as_tibble()

need_cols <- c("deg_overlap", "inference_method", "att_est", "SE",
               "CI_lower", "CI_upper", "ESS_C", "N_T", "SATT", "bias",
               "V_E", "sigma1_extra")
missing <- setdiff(need_cols, names(res))
if (length(missing) > 0)
  stop("Missing columns: ", paste(missing, collapse = ", "))

res2 <- res %>%
  mutate(
    deg_overlap  = factor(deg_overlap, levels = overlap_levels),
    error_satt   = att_est - SATT,
    error_pop    = att_est - ATT_POP_TRUTH,
    covered_satt = (CI_lower <= SATT)          & (SATT           <= CI_upper),
    covered_pop  = (CI_lower <= ATT_POP_TRUTH) & (ATT_POP_TRUTH  <= CI_upper)
  )

summarize_one <- function(df, method_name, sigma1_val) {
  dfm <- df %>%
    filter(inference_method == method_name, sigma1_extra == sigma1_val)

  if (nrow(dfm) == 0) {
    warning("No rows for method=", method_name, " sigma1_extra=", sigma1_val)
    return(tibble())
  }

  dfm %>%
    group_by(deg_overlap) %>%
    summarise(
      N_T         = mean(N_T,    na.rm = TRUE),
      ESS_C       = mean(ESS_C,  na.rm = TRUE),
      avg_SE_hat  = mean(SE,     na.rm = TRUE),
      avg_V_E_hat = mean(V_E,    na.rm = TRUE),
      SE_true_pop  = sd(error_pop,  na.rm = TRUE),
      SE_true_satt = sd(error_satt, na.rm = TRUE),
      Bias_pop     = mean(error_pop,  na.rm = TRUE),
      RMSE_pop     = sqrt(mean(error_pop^2,  na.rm = TRUE)),
      Coverage_pop  = mean(covered_pop,  na.rm = TRUE),
      Coverage_satt = mean(covered_satt, na.rm = TRUE),
      n_runs = sum(!is.na(att_est) & !is.na(SE)),
      .groups = "drop"
    ) %>%
    arrange(deg_overlap) %>%
    mutate(Overlap = overlap_labels[as.character(deg_overlap)]) %>%
    select(Overlap, N_T, ESS_C, avg_SE_hat, avg_V_E_hat,
           SE_true_pop, SE_true_satt,
           Bias_pop, RMSE_pop, Coverage_pop, Coverage_satt, n_runs)
}

sigma1_vals <- c(0, 0.5)
scenario_labels <- c("0" = "Scenario 1 (sigma1=sigma0)", "0.5" = "Scenario 2 (sigma1=sigma0+0.5)")

for (sigma1_val in sigma1_vals) {
  for (m in methods) {
    tab <- summarize_one(res2, m, sigma1_val)
    if (nrow(tab) == 0) next

    scenario_tag <- sprintf("s1extra_%.2f", sigma1_val)
    out_csv <- file.path(out_dir, paste0("table_", m, "_", scenario_tag, ".csv"))
    readr::write_csv(tab, out_csv)

    cat("\n============================================================\n")
    cat("Scenario: sigma1_extra=", sigma1_val, "  (", scenario_labels[as.character(sigma1_val)], ")\n", sep = "")
    cat("Method:   ", m, "\n", sep = "")
    cat("Saved:    ", out_csv, "\n", sep = "")
    cat("------------------------------------------------------------\n")
    print(tab, n = Inf)
  }
}

cat("\nNOTE: Coverage_pop uses population truth ATT=", ATT_POP_TRUTH,
    "; Coverage_satt uses SATT.\n", sep = "")
cat("Expected performance:\n")
cat("  Scenario 1 (sigma1=sigma0): het, alt_common, alt_tt should work; homo should not.\n")
cat("  Scenario 2 (sigma1!=sigma0): only alt_tt should work.\n")
