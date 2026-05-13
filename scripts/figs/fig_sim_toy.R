# Plotting script: scripts/figs/fig_sim_toy.R
# full analysis of simulation results
library(tidyverse)
source("scripts/lib/plot_sim.R")

# res_toy <- readr::read_csv("data/outputs/sims-bias_mse-no-filter-unmatched/toy_combined.csv")
# res_toy <- readr::read_csv("data/outputs/sims-bias_mse-21Nov2025/toy_combined.csv") %>%
#   select(-bal2)
# res_toy_26Apr <- readr::read_csv("data/outputs/sims-bias_mse-26Apr2026/toy_combined.csv") %>%
#   select(-bal2)
res_toy <- readr::read_csv("data/outputs/sims-bias_mse/toy_combined.csv",
                           show_col_types = FALSE)
res_L <- pivot_longer(res_toy,
                      diff:last_col(), names_to="method")

res_L %>% group_by( method ) %>%
  summarise( mean_na = mean( is.na( value ) ),
             SE = sd( value, na.rm=TRUE ),
             sd_tau = sd( true_ATT ) ) %>%
  knitr::kable()




res_toy_L <- prepare_method_comparison_df(res_toy)
res_toy_L

RMSE_plot(res_toy,
          title = "Toy Example",
          xlab = "Value",
          ylab = "Method",
          legend.position = c(0.75, 0.25))

ggsave("figures/sim_toy_results.png", width = 3.5, height = 3.5)


########################################
# diagnostic table view
########################################
res_toy_L <- res_toy %>%
  pivot_longer(
    cols = diff:last_col(), # From diff to the last column
    names_to = "method",
    values_to = "estimate"
  )
res_toy_L

library( simhelpers )
tmp <- res_toy_L %>%
  group_by(method) %>%
  mutate( estimate = estimate - true_ATT,
          true_ATT = 0 ) %>%
  summarise(
    calc_absolute(
      estimates = estimate, true_param = true_ATT,
      criteria = c("bias","stddev", "rmse")
    )
  ) %>%
  arrange(rmse)
tmp

# plot_org_df(org_df_toy,
#             title="Toy Example",
#             xlab="Value",
#             ylab="Method",
#             legend.position=c(0.75, 0.25))
# ggsave("figures/sim_toy_results.png",
#        width=3.5, height=3.5)
