library( CSM )
library(tidyverse)
source("./R/diagnostic_plots.R")
# Test the degree of overlap set by prop_nc_unif
if ( FALSE ){
  sample_dat <-
    gen_one_toy(ctr_dist=0.5,
                prop_nc_unif = 1/3)
  names(sample_dat)
  ggplot( sample_dat, aes( X1, X2, col=as.factor(Z) ) ) +
    geom_point() +
    coord_fixed()
  sample_dat$tau = sample_dat$Y1 - sample_dat$Y0
  skimr::skim( sample_dat )
}


set.seed(123)
# set the same seed to make sure the treated
# is distributed the same as in the simulation
dat_low_overlap <-
  gen_one_toy(ctr_dist = 0.5, prop_nc_unif = 1/3)
plot_low_overlap <-
  create_toy_df_plot(dat_low_overlap) + theme(legend.position = "none")
set.seed(123)
dat_mid_overlap <-
  gen_one_toy(ctr_dist = 0.5, prop_nc_unif = 2/3)
plot_mid_overlap <-
  create_toy_df_plot(dat_mid_overlap) + theme(legend.position = "none")
set.seed(123)
dat_high_overlap <-
  gen_one_toy(ctr_dist = 0.5, prop_nc_unif = 3/3)
plot_high_overlap <-
  create_toy_df_plot(dat_high_overlap)

# Combine plots with a shared legend on the right
combined_plot <- cowplot::plot_grid(
  plot_low_overlap,
  plot_mid_overlap,
  plot_high_overlap,
  ncol = 3, align = 'h',
  rel_widths = c(1, 1, 1.4)
)
print(combined_plot)
ggsave("./figures/sim_toy_3_overlaps.png",
       combined_plot,
       height = 5,
       width = 12)
