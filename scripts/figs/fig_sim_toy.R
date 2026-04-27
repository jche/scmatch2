# Plotting script: scripts/figs/fig_sim_toy.R
# full analysis of simulation results
library(tidyverse)
source("scripts/lib/plot_sim.R")

# res_toy <- readr::read_csv("data/outputs/sims-bias_mse-no-filter-unmatched/toy_combined.csv")
# res_toy <- readr::read_csv("data/outputs/sims-bias_mse/toy_combined.csv")
res_toy <- readr::read_csv("data/outputs/sims-bias_mse-21Nov2025/toy_combined.csv") %>%
  select(-bal2)

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
res_toy_L %>%
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
# plot_org_df(org_df_toy,
#             title="Toy Example",
#             xlab="Value",
#             ylab="Method",
#             legend.position=c(0.75, 0.25))
# ggsave("figures/sim_toy_results.png",
#        width=3.5, height=3.5)
