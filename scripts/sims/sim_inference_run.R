library(tidyverse)
require(mvtnorm)
library( CSM )
source( "scripts/sims/sim_inference_helper.R" )

# Test run
if ( FALSE ){
  toy_naive_low <-
    sim_inference_CSM_A_E(
      dgp_name="toy",
      att0=F,
      R=2,
      prop_nc_unif = 1/3
    )
}


# Run the simulation for inference ----
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
if ( FALSE ){
  source( "scripts/analysis/boot_CSM_simulation_code.R" )
}
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
