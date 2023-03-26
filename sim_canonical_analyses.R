
# full analysis of simulation results


res_kang <- read_csv("sim_canonical_results/kang_spaceship.csv") %>% 
  bind_rows(read_csv("sim_canonical_results/kang_spaceship2.csv")) %>% 
  mutate(sim = "kang")
res_hain <- read_csv("sim_canonical_results/hain_spaceship.csv") %>% 
  mutate(sim = "hain")
res_acic <- read_csv("sim_canonical_results/acic_spaceship.csv") %>% 
  mutate(sim = "acic")

res <- list(res_kang, res_hain, res_acic) %>% 
  map_dfr(function(d) {
    d %>% 
      select(ninf:sim)
  })


# check dropped units -----------------------------------------------------

res %>% 
  ggplot() +
  geom_density(aes(x=ninf)) +
  geom_density(aes(x=ninf_cem), color="red") +
  facet_wrap(~sim, scales="free")



# check rmse stuff --------------------------------------------------------

METHODS <- c("diff", "onenn", "csm_scm", "cem_avg", "bal1", "bal2", 
             "or_lm", "ps_lm",
             "or_bart", "ps_bart",
             "aipw1", "tmle1", 
             "aipw2", "tmle2")

# do ACIC style plots, for all three datasets?
acic_plot <- function(simname, title="", ylab="", xlab="", legend.position="none") {
  res %>% 
    filter(sim == simname) %>% 
    
    pivot_longer(diff:aipw2, names_to="method") %>%
    filter(method %in% METHODS) %>% 
    group_by(method) %>%
    summarize(
      RMSE = sqrt(mean((value-true_ATT)^2)),
      Bias = abs(mean(value-true_ATT)),
    ) %>% 
    mutate(method = fct_reorder(method, RMSE, min)) %>% 
    pivot_longer(c(RMSE, Bias)) %>% 
    
    ggplot(aes(x=method)) +
    geom_point(aes(y=value, pch=name), size=1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
          axis.title.x = element_text(vjust=-3),
          legend.position = legend.position,
          legend.background = element_rect(linetype="solid", linewidth=0.5, color="black")) +
    scale_x_discrete(labels = c(
      "diff" = "diff",
      "tmle1" = "dr-TMLE1",
      "tmle2" = "dr-TMLE2",
      "aipw1" = "dr-AIPW1",
      "aipw2" = "dr-AIPW2",
      "bal1" = "bal-SBW1",
      "bal2" = "bal-SBW2",
      "or_bart" = "or-BART",
      "or_lm" = "or-LM",
      "ps_bart" = "ps-BART",
      "ps_lm" = "ps-LM",
      "onenn" = "match-1NN",
      "cem_avg" = "match-CEM",
      "csm_scm" = expression(bold(match-CSM)))) +
    scale_shape_manual(values = c(1,2)) +
    labs(title = title,
         y = ylab,
         x = xlab,
         pch = "Metric")
}
acic_plot("acic", title="test3", ylab="test", xlab="test2")


require(patchwork)
acic_plot("kang", title="Kang & Schafer", ylab="Value", legend.position=c(0.25, 0.75)) + 
  acic_plot("hain", title="Hainmueller", xlab="Method") + 
  acic_plot("acic", title="ACIC 2016")
  # plot_layout(guides = "collect")
ggsave("writeup/figures/sim_results.png", width=8, height=3.5, units="in")

