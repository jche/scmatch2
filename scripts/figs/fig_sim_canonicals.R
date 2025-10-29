# Plotting script: scripts/figs/fig_sim_canonicals.R
# full analysis of simulation results
library(tidyverse)
source("scripts/analysis/plot_sim.R")

# canonical sims ----------------------------------------------------------
res_kang <- read_csv("data/outputs/sim_canonical_results/kang_spaceship.csv") %>%
  bind_rows(read_csv("data/outputs/sim_canonical_results/kang_spaceship2.csv")) %>%
  mutate(sim = "kang")
res_hain <-
  read_csv("data/outputs/sim_canonical_results/hain_spaceship.csv") %>%
  mutate(sim = "hain")
res_acic <- read_csv("data/outputs/sim_canonical_results/acic_spaceship.csv") %>%
  mutate(sim = "acic")

res <- list(res_kang, res_hain, res_acic) %>%
  map_dfr(function(d) {
    d %>%
      select(ninf:sim)
  })

# check rmse --------------------------------------------------------

# require(patchwork)
acic_plot("kang", title="Kang & Schafer",
          xlab="Method",
          legend.position=c(0.75, 0.25)) +
  acic_plot("hain",ylab="Value", title="Hainmueller") +
  acic_plot("acic", title="ACIC 2016")
ggsave("figures/sim_results.png", width=8, height=3.5, units="in")


# check dropped units -----------------------------------------------------

res %>%
  ggplot() +
  geom_density(aes(x=ninf)) +
  geom_density(aes(x=ninf_cem), color="red") +
  facet_wrap(~sim, scales="free")




