# scripts/analysis/combine_simulation_results.R
# Combine main simulation results with new methods results

library(tidyverse)

# Read main results
main_results <- read_csv("data/outputs/sim_toy_results/toy_comprehensive_run2.csv")

# Read new methods results
new_methods_results <- read_csv("data/outputs/sim_toy_results/toy_new_methods_seq.csv")

# Combine by matching on runid and seed
combined_results <- main_results %>%
  left_join(
    new_methods_results %>%
      select(runid, seed, causal_forest, twang, kbal),
    by = c("runid", "seed")
  )

# Save combined results
write_csv(combined_results, "data/outputs/sim_toy_results/toy_combined_all_methods.csv")

cat("Combined results saved!\n")
cat("Dimensions:", nrow(combined_results), "rows x", ncol(combined_results), "columns\n")

# Quick performance summary
combined_results %>%
  pivot_longer(
    cols = c(diff, bal1, bal2, or_lm, or_bart, ps_lm, ps_bart,
             csm_scm, csm_avg, cem_scm, cem_avg, onenn,
             tmle1, aipw1, tmle3, aipw3,
             causal_forest, twang, kbal),
    names_to = "method",
    values_to = "estimate"
  ) %>%
  group_by(method) %>%
  summarize(
    rmse = sqrt(mean((estimate - true_ATT)^2, na.rm = TRUE)),
    bias = mean(estimate - true_ATT, na.rm = TRUE),
    sd = sd(estimate, na.rm = TRUE),
    n_missing = sum(is.na(estimate))
  ) %>%
  arrange(rmse) %>%
  print(n = Inf)

