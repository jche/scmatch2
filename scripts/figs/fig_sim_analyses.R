# full analysis of simulation results
library(tidyverse)
source("scripts/analysis/plot_sim.R")

# toy sim -----------------------------------------------------------------

res_toy <- readr::read_csv("data/outputs/sim_toy_results/toy_spaceship7.csv")
res_toy <- res_toy %>%  # rename to keep consistency as other results
  rename(tmle2="tmle3",
         aipw2="aipw3")

org_df_toy <-
  res_toy %>%
  pivot_longer(diff:aipw2, names_to="method") %>%
  summarize_bias_rmse()

res_toy_cem_csm <- readr::read_csv("data/outputs/sim_toy_results/toy_spaceship8.csv")

org_df_toy_cem_csm <-
  res_toy_cem_csm %>%
  pivot_longer(csm_scm:cem_avg, names_to="method") %>%
  summarize_bias_rmse()

CEM_CSM_names <- c("csm_scm",
                   "csm_avg",
                   "cem_scm",
                   "cem_avg")

org_df_toy_total <-
  rbind(org_df_toy %>%
        filter(!(method %in%  CEM_CSM_names) ),
      org_df_toy_cem_csm)


plot_org_df(org_df_toy_total,
          title="Toy Example",
          xlab="Value",
          ylab="Method",
          legend.position=c(0.75, 0.25))

# check rmse stuff --------------------------------------------------------
ggsave("figures/sim_toy_results.png",
       width=3.5, height=3.5)


# canonical sims ----------------------------------------------------------
res_kang <- read_csv("sim_canonical_results/kang_spaceship.csv") %>%
  bind_rows(read_csv("sim_canonical_results/kang_spaceship2.csv")) %>%
  mutate(sim = "kang")
res_hain <-
  read_csv("sim_canonical_results/hain_spaceship.csv") %>%
  mutate(sim = "hain")
res_acic <- read_csv("sim_canonical_results/acic_spaceship.csv") %>%
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




