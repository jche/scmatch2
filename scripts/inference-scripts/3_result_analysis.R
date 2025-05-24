# Analysis of 4D Bias Inference Simulation Results
library(tidyverse)
library(knitr)
library(ggplot2)
library(gridExtra)
library(latex2exp)

# Read the combined results file
results <- read_csv("data/outputs/4d_bias_inference_Apr2025/combined_results.csv")

# Function to compute coverage rates
compute_coverage <- results %>%
    mutate(ATT = k * 1.5) %>%
    group_by(deg_overlap) %>%
    summarize(
      n_iterations = n(),
      mean_SATT = mean(SATT, na.rm = TRUE),
      mean_att_est = mean(att_est, na.rm = TRUE),
      mean_bias = mean(bias, na.rm = TRUE),
      mean_SE = mean(SE, na.rm = TRUE),
      true_sd = sd(att_est - bias - ATT, na.rm = TRUE),
      est_sd = mean(SE, na.rm = TRUE),
      true_sd_E = mean(true_sigma * sqrt(1/nt + 1 / ESS_C) , na.rm = TRUE) ,
      true_sd_E_ver_2 = sd(att_est - SATT, na.rm = TRUE),
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
      mean_N_T = mean(N_T, na.rm = TRUE)
    )

# coverage_with_bias_correction <- compute_coverage(results, use_bias_correction = TRUE)
# coverage_without_bias_correction <- compute_coverage(results, use_bias_correction = FALSE)

# 1. Summary table with all metrics
summary_table <- compute_coverage %>%
  select(
    deg_overlap,
    n_iterations,
    mean_SATT,
    mean_att_est,
    mean_bias,
    true_sd,
    est_sd,
    true_sd_E,
    true_sd_E_ver_2,
    est_sd_E,
    mean_V_over_nT,
    true_var,
    coverage_rate,
    coverage_rate_wo_bias_correction,
    avg_CI_width,
    mean_ESS_C,
    mean_N_T
  )

print(kable(summary_table,
            caption = "Summary of 4D Bias Inference Simulation Results (with bias correction)",
            digits = 3))

# 2. Compare coverage with and without bias correction
coverage_comparison <- bind_rows(
  coverage_with_bias_correction %>% mutate(bias_correction = "With Bias Correction"),
  coverage_without_bias_correction %>% mutate(bias_correction = "Without Bias Correction")
) %>%
  select(deg_overlap, bias_correction, coverage_rate, avg_CI_width, mean_bias)

# 3. Analysis for checking theoretical assumptions

# 3.1 Calculate SATT based on dimensions (k=4 in this case)
# Based on paper, in k dimensions, SATT should be 1.5 * k
k <- 4
theoretical_SATT <- 1.5 * k # = 6 for k=4

# Verify the average SATT matches our theoretical value
actual_SATT <- mean(results$SATT, na.rm = TRUE)
cat("Theoretical SATT (1.5 * k =", theoretical_SATT, ") vs. Actual mean SATT:", actual_SATT, "\n")

# 3.2 Check for V and V_E consistency

# Calculate true V_E using variance of estimation error
V_E_empirical <- results %>%
  group_by(deg_overlap) %>%
  summarize(
    # Method 1: Var(att_est - bias)
    V_E_method1 = var(att_est - bias, na.rm = TRUE) * mean(N_T, na.rm = TRUE),
    # Method 2: sigma^2 * (1/n_T + 1/ESS)
    # Assuming sigma^2 is approximately the pooled variance of residuals
    mean_N_T = mean(N_T, na.rm = TRUE),
    mean_ESS_C = mean(ESS_C, na.rm = TRUE),
    # Calculate estimated V_E from SE^2 * n_T since SE = sqrt(V/n_T)
    estimated_V_E = mean(SE^2, na.rm = TRUE) * mean(N_T, na.rm = TRUE),
    # Ratio of estimated to empirical
    ratio_method1 = estimated_V_E / V_E_method1
  )

# Checking bias significance by dimension
if (k <= 2) {
  # For k ≤ 2, bias should be negligible
  cat("Since k =", k, "≤ 2, bias should be negligible.\n")
  cat("Mean bias across all iterations:", mean(results$bias, na.rm = TRUE), "\n")
} else {
  # For k > 2, bias is significant and needs correction
  cat("Since k =", k, "> 2, bias is significant and needs correction.\n")
  cat("Mean bias across all iterations:", mean(results$bias, na.rm = TRUE), "\n")
}

# 4. Plotting functions

# 4.1 Plot coverage rates with and without bias correction
plot_coverage_comparison <- function() {
  ggplot(coverage_comparison, aes(x = deg_overlap, y = coverage_rate, fill = bias_correction)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    ylim(0, 1) +
    labs(
      title = "Coverage Rates by Degree of Overlap",
      subtitle = "Comparison with and without bias correction",
      x = "Degree of Overlap",
      y = "Coverage Rate",
      fill = "Method"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# 4.2 Plot CI width with and without bias correction
plot_CI_width_comparison <- function() {
  ggplot(coverage_comparison, aes(x = deg_overlap, y = avg_CI_width, fill = bias_correction)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      title = "Average CI Width by Degree of Overlap",
      subtitle = "Comparison with and without bias correction",
      x = "Degree of Overlap",
      y = "Average CI Width",
      fill = "Method"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# 4.3 Plot the distribution of normalized estimation errors
plot_normalized_errors <- function() {
  # Calculate normalized errors: (att_est - bias - SATT) / SE
  results_with_normalized <- results %>%
    mutate(
      normalized_error = (att_est - bias - SATT) / SE
    )

  # Histogram with normal density overlay
  ggplot(results_with_normalized, aes(x = normalized_error)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = "red", size = 1) +
    labs(
      title = "Distribution of Normalized Estimation Errors",
      subtitle = "With theoretical N(0,1) density overlay",
      x = "Normalized Error",
      y = "Density"
    ) +
    theme_minimal()
}

# 4.4 Plot bias by degree of overlap
plot_bias_by_overlap <- function() {
  results %>%
    group_by(deg_overlap) %>%
    summarize(
      mean_bias = mean(bias, na.rm = TRUE),
      se_bias = sd(bias, na.rm = TRUE) / sqrt(n())
    ) %>%
    ggplot(aes(x = deg_overlap, y = mean_bias)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_errorbar(aes(ymin = mean_bias - 1.96 * se_bias, ymax = mean_bias + 1.96 * se_bias), width = 0.2) +
    labs(
      title = "Average Bias by Degree of Overlap",
      x = "Degree of Overlap",
      y = "Average Bias"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# 4.5 QQ plot for normality check
plot_qq <- function() {
  results_with_normalized <- results %>%
    mutate(
      normalized_error = (att_est - bias - SATT) / SE
    )

  ggplot(results_with_normalized, aes(sample = normalized_error)) +
    stat_qq() +
    stat_qq_line(color = "red") +
    labs(
      title = "QQ Plot of Normalized Estimation Errors",
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_minimal()
}

# 4.6 Plot estimated vs empirical standard errors
plot_se_comparison <- function() {
  se_comparison <- bind_rows(
    coverage_with_bias_correction %>%
      select(deg_overlap, sd_value = theoretical_sd) %>%
      mutate(type = "Theoretical SE"),
    coverage_with_bias_correction %>%
      select(deg_overlap, sd_value = emp_sd) %>%
      mutate(type = "Empirical SE")
  )

  ggplot(se_comparison, aes(x = deg_overlap, y = sd_value, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      title = "Comparison of Theoretical vs. Empirical Standard Errors",
      x = "Degree of Overlap",
      y = "Standard Error",
      fill = "Type"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Generate all plots
p1 <- plot_coverage_comparison()
p2 <- plot_CI_width_comparison()
p3 <- plot_normalized_errors()
p4 <- plot_bias_by_overlap()
p5 <- plot_qq()
p6 <- plot_se_comparison()

# Print plots
p1
p2
p3
p4
p5
p6

# Create a comparison for V_E consistency
V_E_comparison <- V_E_empirical %>%
  select(deg_overlap, V_E_method1, estimated_V_E, ratio_method1)

print(kable(V_E_comparison,
            caption = "Comparison of Empirical vs. Estimated V_E",
            digits = 3))

# Save results to files
write.csv(summary_table, "data/outputs/4d_bias_inference_Apr2025/summary_table.csv", row.names = FALSE)
write.csv(coverage_comparison, "data/outputs/4d_bias_inference_Apr2025/coverage_comparison.csv", row.names = FALSE)
write.csv(V_E_comparison, "data/outputs/4d_bias_inference_Apr2025/V_E_comparison.csv", row.names = FALSE)

# Additional code for LaTeX output
latex_summary_table <- kable(summary_table, format = "latex",
                             caption = "Summary of 4D Bias Inference Simulation Results (with bias correction)",
                             digits = 3)

cat("LaTeX Summary Table:\n")
cat(latex_summary_table)

# Save a combined plot for the paper
combined_plot <- grid.arrange(p1, p3, p5, p6, ncol = 2)
ggsave("data/outputs/4d_bias_inference_Apr2025/combined_plots.pdf", combined_plot, width = 12, height = 10)
