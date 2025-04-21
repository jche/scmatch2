
library(tidyverse)
require(mvtnorm)
library( CSM )
source( here::here( "scripts/inference-scripts/0_sim_inference_utils.R" ) )

# Test run ------
if ( FALSE ){
  R = 2
  toy_naive_low_A_E <-
    sim_inference_CSM_A_E(
      R=R,
      prop_nc_unif = 1/3,
      seed = c(123 + 1:R * 2)
    )
  # att_true is the SATT without noise
  # att_est is the SATT plus noise

  toy_naive_low_A_E

  toy_naive_low_OR <-
    sim_inference_CSM_OR(
      R=R,
      prop_nc_unif = 1/3,
      seed = c(123 + 1:R * 2)
    )
  toy_naive_low_OR
}

# New code -----
#' Generalized function to run simulations for sim_inference_CSM_A_E or sim_inference_CSM_OR
#' Saves results to the specified filename
run_sim_inference <- function(R = 10, method = c("A_E", "OR"), FNAME,
                              toy_ctr_dist = 0.5,
                              scaling = 8,
                              true_sigma = 0.5,
                              nc = 500,
                              parallel = TRUE) {
  method <- match.arg(method)

  ### 5 degrees of freedom
  prop_nc_unif_values <- seq(0.2, 1, by = 0.2)
  deg_overlap_labels <- c("very_low", "low", "mid", "high", "very_high")

  # Select appropriate function
  sim_function <- if (method == "A_E") sim_inference_CSM_A_E else sim_inference_CSM_OR

  for (i in seq_along(prop_nc_unif_values)) {
    cat("Simulation", i, "(", method, ")\n")
    sim_result <- sim_function(
      R = R,
      prop_nc_unif = prop_nc_unif_values[i],
      scaling = scaling,
      toy_ctr_dist = toy_ctr_dist,
      true_sigma = true_sigma,
      nc = nc,
      seed = c(123 + i * 2),
      parallel = parallel
    )

    sim_result$deg_overlap <- deg_overlap_labels[i]

    save_res_to_csv(sim_result, FNAME = FNAME)
  }

  cat("Simulation complete for method:", method, "\n")

  invisible(0)
}

# Run the simulation ----

if (TRUE) {
  R = 500

  # for (method in c("A_E", "OR")) {
  for (method in c("OR")) {
    file_folder <- here::here(paste0("data/outputs/", method, "-overlap-by-prop-unif"))
    dir.create(file_folder, showWarnings = FALSE, recursive = TRUE)
    FNAME = here::here(
      paste0("data/outputs/", method, "-overlap-by-prop-unif/",
             method, "_toy_low_mid_high_R=", R, ".csv")
    )

    file.remove(FNAME)

    # Run simulation
    tictoc::tic()
    run_sim_inference(R = R,
                      method = method,
                      FNAME = FNAME,
                      scaling = 8,
                      true_sigma = 0.5,
                      nc = 500)
    tictoc::toc()

    # Check saved file
    rs = read_csv(FNAME)
    skimr::skim(rs)

    cat("Results saved to", FNAME, "\n")
  }
}
