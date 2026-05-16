
# Run just a few scenarios locally to see if coverage is good.


suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(furrr)
  library(tictoc)
  library(here)
  library(tidyverse)

})

N_ITER      <- 100
OUTPUT_NAME <- "sims-quick"
SAVE_OUTPUT <- TRUE

# ── Load utilities ────────────────────────────────────────────────────────────
source(here::here("scripts/sims-variance-multi/0_utils.R"))

save_path = NULL
if (SAVE_OUTPUT) {
  paths <- get_sim_paths(OUTPUT_NAME)
  message(sprintf("Saving output to: %s", paths$individual_dir))
  save_path = paths$individual_dir
}



DESIGN_GRID <- expand_grid(
  overlap_label = "mid",
  error_type    = c("homo", "het"),
  sigma1_extra  = c(0, 2),
  tx_type       = c("het", "constant"),
  k_match = c( 2 )
) %>%
  mutate(
    common_label = if_else(sigma1_extra == 0, "common", "no_common"),
    # factor ordering for plots
    overlap_label = factor(overlap_label, levels = OVERLAP_LABELS),
    error_type    = factor(error_type,    levels = c("homo", "het")),
    common_label  = factor(common_label,  levels = c("common", "no_common"))
  )
DESIGN_GRID


# ── Run in parallel ───────────────────────────────────────────────────────────
n_workers <- parallel::detectCores() - 2
message(sprintf(
  "Running %d iteration(s) × %d design cells in parallel (%d workers)  [save=%s]",
  N_ITER, nrow(DESIGN_GRID), n_workers, SAVE_OUTPUT
))

plan(multisession, workers = n_workers)
tic()

all_res <- future_map( 1:N_ITER,
                       .f = sim_master_multi,
                       caliper       = 0.1,
                       save_path = save_path,
                       overwrite = TRUE,
                       .options = furrr_options(seed = NULL),
                       .progress = TRUE
)

toc()
plan(sequential)


if ( FALSE ) {
  #Testing
  sim_master_multi( 4, save_path=save_path, overwrite = FALSE )
  debugonce( sim_master_multi )
  sim_master_multi( 11, save_path=save_path, overwrite = FALSE )

}


cat("══════════════════════════════════════════════════════════════════════\n")
cat( "\n\nQuick Simulation complete\n" )
cat("══════════════════════════════════════════════════════════════════════\n")

res <- collect_sim_results( OUTPUT_NAME )

res

cov <- res %>%
  group_by( error_type, sigma1_extra, common_label, k_match, tx_type, method ) %>%
  mutate( cover = SATT <= att_est+2*SE & SATT >= att_est-2*SE ) %>%
  summarize( R = n(),
             mean_S0 = sqrt( mean( S0_sq ) ),
             cover = mean( cover ), .groups = "drop" )

cov %>%
  arrange( common_label, method )


ggplot( cov, aes( error_type, cover, color=method, group=method ) ) +
  facet_grid(  tx_type ~ common_label ) +
  geom_point() +
  geom_line() +
  geom_hline( yintercept = 0.95 )



res %>%
  pivot_longer( cols = c("S0_sq", "S1_sq") ) %>%
  mutate( value = sqrt( value ) ) %>%
  ggplot( aes( error_type, value, col = name ) ) +
  facet_grid(   ~ common_label ) +
  geom_boxplot()




