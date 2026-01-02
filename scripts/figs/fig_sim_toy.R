# Plotting script: scripts/figs/fig_sim_toy.R
# full analysis of simulation results
library(tidyverse)
source("scripts/analysis/plot_sim.R")

res_toy <- readr::read_csv("data/outputs/sims-bias_mse/toy_combined.csv")

org_df_toy <- res_toy %>%
  select(-bal2) %>%
  pivot_longer(
    cols = diff:last_col(),  # From diff to the last column
    names_to = "method"
  ) %>%
  summarize_bias_rmse()

plot_org_df(org_df_toy,
            title="Toy Example",
            xlab="Value",
            ylab="Method",
            legend.position=c(0.75, 0.25))

ggsave("figures/sim_toy_results.png",
       width=3.5, height=3.5)

# # toy sim on sim_toy with all methods  -----------------------------------------------------------------
# res_toy <- readr::read_csv("data/outputs/sim_toy_results/toy_combined_all_methods.csv")
# res_toy <- res_toy %>%  # rename to keep consistency as other results
#   rename(tmle2="tmle3",
#          aipw2="aipw3")
#
# org_df_toy <- res_toy %>%
#   pivot_longer(
#     cols = diff:last_col(),  # From diff to the last column
#     names_to = "method"
#   ) %>%
#   summarize_bias_rmse()
#
# plot_org_df(org_df_toy,
#             title="Toy Example",
#             xlab="Value",
#             ylab="Method",
#             legend.position=c(0.75, 0.25))
#
# ggsave("figures/sim_toy_results_seeded_all_method.png",
#        width=3.5, height=3.5)


# # toy sim -----------------------------------------------------------------
#
# res_toy <- readr::read_csv("data/outputs/sim_toy_results/toy_spaceship7.csv")
# res_toy <- res_toy %>%  # rename to keep consistency as other results
#   rename(tmle2="tmle3",
#          aipw2="aipw3")
#
# org_df_toy <-
#   res_toy %>%
#   pivot_longer(diff:aipw2, names_to="method") %>%
#   summarize_bias_rmse()
#
# res_toy_cem_csm <- readr::read_csv("data/outputs/sim_toy_results/toy_spaceship8.csv")
#
# org_df_toy_cem_csm <-
#   res_toy_cem_csm %>%
#   pivot_longer(csm_scm:cem_avg, names_to="method") %>%
#   summarize_bias_rmse()
#
# CEM_CSM_names <- c("csm_scm",
#                    "csm_avg",
#                    "cem_scm",
#                    "cem_avg")
#
# org_df_toy_total <-
#   rbind(org_df_toy %>%
#         filter(!(method %in%  CEM_CSM_names) ),
#       org_df_toy_cem_csm)
#
#
# plot_org_df(org_df_toy,
#           title="Toy Example",
#           xlab="Value",
#           ylab="Method",
#           legend.position=c(0.75, 0.25))
#
# # check rmse stuff --------------------------------------------------------
# ggsave("figures/sim_toy_results.png",
#        width=3.5, height=3.5)
