
# Run just a few scenarios locally to see if coverage is good.


suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(furrr)
  library(tictoc)
  library(here)
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
                       save_path = save_path,
                       overwrite = FALSE,
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


cat( "Quick Simulation complete\n" )

all_res = all_res[ !is.na( all_res ) ]
res = bind_rows( all_res )
res

cov <- res %>%
  group_by( error_type, sigma1_extra, common_label, k_match,method ) %>%
  mutate( cover = SATT <= att_est+2*SE & SATT >= att_est-2*SE ) %>%
  summarize( cover = mean( cover ), .groups = "drop" )



