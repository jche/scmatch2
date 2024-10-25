
library(tidyverse)
require(mvtnorm)
library( CSM )
source( here::here( "scripts/sims/sim_inference_helper.R" ) )

# Test run
if ( FALSE ){
  R = 2
  toy_naive_low <-
    sim_inference_CSM_A_E(
      dgp_name="toy",
      att0=F,
      R=R,
      prop_nc_unif = 1/3,
      seed = c(123 + 1:R * 2)
    )
  # att_true is the SATT without noise
  # att_est is the SATT plus noise
}


run_sim_inference_A_E <- function(R = 10, FNAME) {
  ### 3 degrees of freedoms
  # prop_nc_unif_values <- c(1/3, 2/3, 3/3)
  # deg_overlap_labels <- c("low", "mid", "high")

  ### 5 degrees of freedoms
  prop_nc_unif_values <- seq(0.2, 1, by = 0.2)
  deg_overlap_labels <- c("very_low", "low", "mid", "high", "very_high")

  for (i in seq_along(prop_nc_unif_values)) {
    cat( "Simulation", i, "\n" )
    sim_result <- sim_inference_CSM_A_E(
      dgp_name = "toy",
      att0 = FALSE,
      R = R,
      prop_nc_unif = prop_nc_unif_values[i],
      seed = c(123 + i*2),
      parallel = TRUE
    )

    sim_result$deg_overlap <- deg_overlap_labels[i]

    save_res_to_csv(sim_result, FNAME = FNAME)
  }

  cat( "Simulation complete\n" )
}



if ( TRUE ) {
  tictoc::tic()
  R = 500
  FNAME =
    here::here(
      paste0("data/outputs/A-E-overlap-by-prop-unif/",
             "A_E_toy_low_mid_high_R=",R,".csv")
    )
  file.remove( FNAME )
  run_sim_inference_A_E(R=R, FNAME = FNAME)
  tictoc::toc()

  rs = read_csv( FNAME )
  skimr::skim( rs )

  cat("Results saved to", FNAME, "\n" )
}
