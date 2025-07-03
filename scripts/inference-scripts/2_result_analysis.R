library(tidyverse)
library(knitr)
library(ggplot2)
library(gridExtra)
library(latex2exp)

# Combine the results
source(here::here("scripts/inference-scripts/0_sim_inference_utils.R"))
collect_results_to_csv()

# Read the combined results file
paths <- get_sim_paths()
combined_csv_file <- paths$combined_csv
results <- read_csv(combined_csv_file)

tmp <- results %>%
  mutate(ATT = k * 1.5) %>%
  filter(deg_overlap == "high", inference_method=="pooled")


# Function to compute coverage rates
diagnostic_table <- results %>%
    mutate(ATT = k * 1.5) %>%
    group_by(deg_overlap, inference_method) %>%
    summarize(
      n_iterations = n(),
      mean_SATT = mean(SATT, na.rm = TRUE),
      mean_att_est = mean(att_est, na.rm = TRUE),
      mean_bias = mean(bias, na.rm = TRUE),
      mean_SE = mean(SE, na.rm = TRUE),
      true_sd = sd(att_est - bias - ATT, na.rm = TRUE),
      est_sd = mean(SE, na.rm = TRUE),
      true_sd_E_ver_2 = mean(true_sigma * sqrt(1/nt + 1 / ESS_C) , na.rm = TRUE) ,
      true_sd_E = sd(att_est - SATT, na.rm = TRUE),
      est_sd_E = mean(sigma_hat * sqrt(1/nt + 1 / ESS_C) , na.rm = TRUE) ,
      # est_sd_E_ver_2 = mean(sqrt(V_E), na.rm = TRUE) , # Exactly same as est_sd_E
      coverage_rate = mean((CI_lower - bias <= ATT) & (ATT <= CI_upper - bias), na.rm = TRUE),
      coverage_rate_wo_bias_correction = mean((CI_lower <= ATT) & (ATT <= CI_upper), na.rm = TRUE),
      avg_CI_width = mean((CI_upper - CI_lower), na.rm = TRUE),
      # Average estimated V/n_T (from SE^2)
      mean_V_over_nT = mean(SE^2, na.rm = TRUE),
      # Empirical variance
      true_var = var(att_est - bias - SATT, na.rm = TRUE) * 10,
      # Mean effective sample size of controls
      mean_ESS_C = mean(ESS_C, na.rm = TRUE),
      # Number of treated units
      mean_N_T = mean(N_T, na.rm = TRUE),
      # Mean overlap statistics
      mean_avg_shared_controls = mean(avg_shared_controls, na.rm = TRUE),
      mean_p75_shared_controls = mean(p75_shared_controls, na.rm = TRUE),
      mean_avg_shared_treated = mean(avg_shared_treated, na.rm = TRUE),
      mean_p75_shared_treated = mean(p75_shared_treated, na.rm = TRUE),
      .groups = 'drop'
    )

# Main LaTeX table for paper
main_table <- diagnostic_table %>%
  mutate(control_reuse = round(mean_avg_shared_controls, 1)) %>%
  select(deg_overlap,
         inference_method,
         coverage_rate,
         est_sd,
         control_reuse
         )

print(kable(main_table,
            caption = "Che et al. DGP: Coverage Performance Across Overlap Scenarios",
            digits = 3))


# 2. Table to check for V and V_E consistency

# Calculate true V_E using variance of estimation error
supporting_table <- diagnostic_table %>%
  select(
    deg_overlap,
    inference_method,
    est_sd,
    true_sd_E,
    est_sd_E,
    coverage_rate,
    coverage_rate_wo_bias_correction,
    avg_CI_width,
    mean_ESS_C,
  )

print(kable(supporting_table,
            caption = "Che et al. DGP: Additional Statistics",
            digits = 3))



#
# # 1. Summary table with all metrics
# summary_table <- results_summarized %>%
#   dplyr::select(
#     deg_overlap,
#     inference_method,
#     n_iterations,
#     mean_SATT,
#     mean_att_est,
#     mean_bias,
#     true_sd,
#     est_sd,
#     true_sd_E,
#     est_sd_E,
#     coverage_rate,
#     coverage_rate_wo_bias_correction,
#     avg_CI_width,
#     mean_ESS_C,
#     mean_N_T,
#     mean_avg_shared_controls,
#     mean_p75_shared_controls,
#     mean_avg_shared_treated,
#     mean_p75_shared_treated
#   )
#
# print(kable(summary_table,
#             caption = "Summary of 4D Bias Inference Simulation Results (with bias correction)",
#             digits = 3))
