library(tidyverse)
require(mvtnorm)

library( CSM )
source( "scripts/analysis/boot_CSM_simulation_code.R" )

save_res_to_csv<-
  function(curr_res,
           FNAME){
    if (file.exists(FNAME)) {
      write_csv(curr_res, FNAME, append=TRUE)
    } else {
      write_csv(curr_res, FNAME)
    }
  } # save_res_to_csv




sample_dat <-
  gen_one_toy(ctr_dist=0.5,
              prop_nc_unif = 1/3)
names(sample_dat)
ggplot( sample_dat, aes( X1, X2, col=as.factor(Z) ) ) +
  geom_point() +
  coord_fixed()
sample_dat$tau = sample_dat$Y1 - sample_dat$Y0
skimr::skim( sample_dat )

## plot the
source("./R/diagnostic_plots.R")
# Generate data and create individual plots
dat_low_overlap <-
  gen_one_toy(ctr_dist = 0.5, prop_nc_unif = 1/3)
plot_low_overlap <-
  create_toy_df_plot(dat_low_overlap) + theme(legend.position = "none")

dat_mid_overlap <-
  gen_one_toy(ctr_dist = 0.5, prop_nc_unif = 2/3)
plot_mid_overlap <-
  create_toy_df_plot(dat_mid_overlap) + theme(legend.position = "none")

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
  rel_widths = c(1, 1, 1.4)  # Adjust widths if needed
)
print(combined_plot)
ggsave("./writeup/figures/sim_toy_3_overlaps.png",
       combined_plot,
       height = 5,
       width = 12)

# Test simulation driver ----

if ( FALSE ) {
  tst <-
    boot_CSM(dgp_name="toy",
             att0=F,
             I=5,
             B=10,
             mu_model="linear",
             boot_mtd="A-E",
             n_split=1,
             kang_true = FALSE,
             toy_ctr_dist = 0.5)

}

# Run the simulation for inference ----
## With the sim_inference_CSM_A_E code, where
#   the inference is run on
run_sim_inference_A_E <- function(R = 10){

  # R <- 2; source( "scripts/analysis/boot_CSM_simulation_code.R" )
  toy_naive_low <-
    sim_inference_CSM_A_E(
      dgp_name="toy",
      att0=F,
      R=R,
      # toy_ctr_dist=0.5,
      prop_nc_unif = 1/3
    )
  toy_naive_low$deg_overlap <- "low"
  save_res_to_csv(toy_naive_low,
                  FNAME = FNAME)


  toy_naive_mid<-
    sim_inference_CSM_A_E(
      dgp_name="toy",
      att0=F,
      R=R,
      # toy_ctr_dist=0.3,
      prop_nc_unif = 2/3
    )
  toy_naive_mid$deg_overlap <- "mid"
  save_res_to_csv(toy_naive_mid,
                  FNAME = FNAME)


  toy_naive_high<-
    sim_inference_CSM_A_E(
      dgp_name="toy",
      att0=F,
      R=R,
      # toy_ctr_dist=0.1,
      prop_nc_unif = 3/3
    )
  toy_naive_high$deg_overlap <- "high"
  save_res_to_csv(toy_naive_high,
                  FNAME = FNAME)
}
# source( "scripts/analysis/boot_CSM_simulation_code.R" )
{
  tictoc::tic()
  R = 1000
  FNAME =
    here::here(
      paste0("data/outputs/A-E-overlap-by-prop-unif/",
             "A_E_toy_low_mid_high_R=",R,".csv")
    )
  run_sim_inference_A_E(R=R)
  tictoc::toc()
}
