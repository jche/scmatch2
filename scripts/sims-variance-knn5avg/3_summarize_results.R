# scripts/sims-variance/3_summarize_results.R
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(here)
  library(stringr)
})

# input
in_csv <- here::here("data/outputs/sims-variance-knn5avg/combined_results.csv")


# output dir
out_dir <- here::here("tables/sims-variance-knn5avg")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# methods to report
methods <- c("pooled", "pooled_het", "pooled_V_E_CI", "pooled_het_V_E_CI")

# overlap order + pretty names
overlap_levels <- c("very_low", "low", "mid", "high", "very_high")
overlap_labels <- c(
  very_low  = "Very Low",
  low       = "Low",
  mid       = "Medium",
  high      = "High",
  very_high = "Very High"
)

# population truth requested by you
ATT_POP_TRUTH <- 3

stopifnot(file.exists(in_csv))
res <- readr::read_csv(in_csv, show_col_types = FALSE) %>% as_tibble()

# quick sanity checks
need_cols <- c("deg_overlap", "inference_method", "att_est", "SE", "CI_lower", "CI_upper",
               "ESS_C", "N_T", "SATT", "bias")
missing <- setdiff(need_cols, names(res))
if (length(missing) > 0) {
  stop("Missing required columns in combined results: ", paste(missing, collapse = ", "))
}

# derived columns
res2 <- res %>%
  mutate(
    deg_overlap = factor(deg_overlap, levels = overlap_levels),

    # --- errors vs sample truth (SATT) ---
    error_satt = att_est - SATT,
    noise_satt = error_satt - bias,

    # --- errors vs population truth (3) ---
    error_pop = att_est - ATT_POP_TRUTH,
    noise_pop = error_pop - bias,

    # coverage vs population truth
    covered_pop = (CI_lower <= ATT_POP_TRUTH) & (ATT_POP_TRUTH<= CI_upper),

    # coverage vs SATT
    covered_satt = (CI_lower <= SATT) & (SATT <= CI_upper)
  )

summarize_one_method <- function(df, method_name) {
  dfm <- df %>% filter(inference_method == method_name)

  if (nrow(dfm) == 0) {
    warning("No rows for method: ", method_name)
    return(tibble())
  }

  # ---- method-specific rules you requested ----
  if (method_name %in% c("pooled", "pooled_het")) {
    cov_col <- "covered_pop"
    # se_true_fun <- function(d) sd(d$att_est, na.rm = TRUE)
    se_true_fun <- function(d) sd(d$error_pop , na.rm = TRUE)
    bias_col <- "error_pop"
  } else if (method_name %in% c("pooled_V_E_CI", "pooled_het_V_E_CI")) {
    cov_col <- "covered_satt"
    se_true_fun <- function(d) sd(d$error_satt, na.rm = TRUE)
    bias_col <- "error_satt"
  } else {
    stop("Unknown method_name: ", method_name)
  }

  tab <- dfm %>%
    group_by(deg_overlap) %>%
    summarise(
      N_T = mean(N_T, na.rm = TRUE),
      N_C_tilde = mean(ESS_C, na.rm = TRUE),

      avg_SE_hat = mean(SE, na.rm = TRUE),

      # method-specific "true SE"
      avg_SE_true = se_true_fun(cur_data()),

      # Bias/RMSE aligned with the same truth used above
      Bias = mean(.data[[bias_col]], na.rm = TRUE),
      RMSE = sqrt(mean((.data[[bias_col]])^2, na.rm = TRUE)),

      Coverage = mean(.data[[cov_col]], na.rm = TRUE),

      n_runs = sum(!is.na(att_est) & !is.na(SE)),
      .groups = "drop"
    ) %>%
    arrange(deg_overlap) %>%
    mutate(
      Overlap = overlap_labels[as.character(deg_overlap)]
    ) %>%
    select(
      Overlap,
      N_C_tilde,
      avg_SE_hat,
      avg_SE_true,
      Bias,
      RMSE,
      Coverage,
      n_runs
    )

  tab
}

for (m in methods) {
  tab_m <- summarize_one_method(res2, m)
  if (nrow(tab_m) == 0) next

  out_csv <- file.path(out_dir, paste0("table_knn5avg_", m, ".csv"))

  readr::write_csv(tab_m, out_csv)

  cat("\n============================================================\n")
  cat("Method: ", m, "\n", sep = "")
  cat("Saved:  ", out_csv, "\n", sep = "")
  cat("------------------------------------------------------------\n")
  print(tab_m, n = Inf)
}

cat("\nNOTE:\n")
cat(" - pooled / pooled_het: Coverage uses population truth ATT = ", ATT_POP_TRUTH, " (covered_pop)\n", sep = "")
cat(" - pooled_V_E_CI / pooled_het_V_E_CI: Coverage uses SATT (covered_satt)\n", sep = "")
