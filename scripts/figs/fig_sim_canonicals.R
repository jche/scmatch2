# Plotting script: scripts/figs/fig_sim_canonicals.R
# full analysis of simulation results
library(tidyverse)
source("scripts/analysis/plot_sim.R")

# canonical sims-bias_mse ----------------------------------------------------------
res_acic <- read_csv("data/outputs/sims-bias_mse/acic_combined.csv")
res_hain <- read_csv("data/outputs/sims-bias_mse/hainmueller_combined.csv")
res_kang <- read_csv("data/outputs/sims-bias_mse/kang_combined.csv")

res <- list(res_kang, res_hain, res_acic) %>%
  map_dfr(function(d) {
    # d %>%
    #   select(ninf:last_col())
    d
  })

# check rmse --------------------------------------------------------

require(patchwork)
acic_plot_sim_type("kang", title="Kang & Schafer",
          xlab="Method",
          legend.position=c(0.75, 0.25)) +
  acic_plot_sim_type("hainmueller",ylab="Value", title="Hainmueller") +
  acic_plot_sim_type("acic", title="ACIC 2016")
ggsave("figures/sim_canonical_results.png", width=8, height=3.5, units="in")


# # check dropped units -----------------------------------------------------
#
# res %>%
#   ggplot() +
#   geom_density(aes(x=ninf)) +
#   geom_density(aes(x=ninf_cem), color="red") +
#   facet_wrap(~sim, scales="free")
#



