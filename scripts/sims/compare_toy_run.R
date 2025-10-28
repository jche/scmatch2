library(tidyverse)
source("scripts/analysis/plot_sim.R")

# Load both datasets
comp_run <- read_csv("data/outputs/sim_toy_results/toy_comprehensive_run.csv")
spaceship7 <- read_csv("data/outputs/sim_toy_results/toy_spaceship7.csv")

# Check dimensions
cat("Dimensions:\n")
cat("comp_run:", dim(comp_run), "\n")
cat("spaceship7:", dim(spaceship7), "\n\n")

# Check column names differences
cat("Columns in comp_run but not spaceship7:\n")
print(setdiff(names(comp_run), names(spaceship7)))
cat("\nColumns in spaceship7 but not comp_run:\n")
print(setdiff(names(spaceship7), names(comp_run)))

# Rename tmle3/aipw3 to tmle2/aipw2 for both
comp_run <- comp_run %>%
  rename(tmle2 = "tmle3", aipw2 = "aipw3")

spaceship7 <- spaceship7 %>%
  rename(tmle2 = "tmle3", aipw2 = "aipw3")

# Create summary dataframes
org_df_comp <- comp_run %>%
  pivot_longer(diff:aipw2, names_to = "method") %>%
  summarize_bias_rmse()

org_df_spaceship7 <- spaceship7 %>%
  pivot_longer(diff:aipw2, names_to = "method") %>%
  summarize_bias_rmse()

# Join by method and name to compare exact values
comparison <- org_df_comp %>%
  rename(value_comp = value) %>%
  full_join(
    org_df_spaceship7 %>% rename(value_spaceship7 = value),
    by = c("method", "name")
  ) %>%
  mutate(
    difference = value_comp - value_spaceship7,
    abs_difference = abs(difference),
    pct_difference = 100 * difference / value_spaceship7
  ) %>%
  arrange(desc(abs_difference))

# Print comparison
cat("\n=== Comparison of comp_run vs spaceship7 ===\n")
print(comparison, n = Inf)

# Summary statistics
cat("\n=== Summary of Differences ===\n")
cat("Mean absolute difference:", mean(comparison$abs_difference, na.rm = TRUE), "\n")
cat("Max absolute difference:", max(comparison$abs_difference, na.rm = TRUE), "\n")
cat("Correlation between values:",
    cor(comparison$value_comp, comparison$value_spaceship7, use = "complete.obs"), "\n")

# Check for any major discrepancies (>10% difference)
major_diffs <- comparison %>%
  filter(abs(pct_difference) > 10 | is.na(pct_difference) | is.infinite(pct_difference))

if (nrow(major_diffs) > 0) {
  cat("\n=== Methods with >10% difference ===\n")
  print(major_diffs)
} else {
  cat("\nNo methods with >10% difference - datasets appear consistent!\n")
}

# Plot both for visual comparison
plot_comp <- plot_org_df(org_df_comp,
                         title = "Comprehensive Run",
                         xlab = "Value",
                         ylab = "Method",
                         legend.position = c(0.75, 0.25))

plot_spaceship7 <- plot_org_df(org_df_spaceship7,
                               title = "Spaceship7",
                               xlab = "Value",
                               ylab = "Method",
                               legend.position = c(0.75, 0.25))

# Side by side comparison (if you have patchwork or gridExtra)
library(patchwork)
plot_comp + plot_spaceship7
