# scripts/sims-variance-multi/2_collect.R
#
# Collect per-iteration .rds files into a single combined CSV + RDS.
# Logic lives in collect_sim_results() in 0_utils.R.
#
# Usage:
#   Rscript 2_collect.R [output_name]
#   output_name  default: sims-variance-multi

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(here)
  library( tidyverse )
})

args        <- commandArgs(trailingOnly = TRUE)
output_name <- if (length(args) >= 1) args[[1]] else "sims-variance-multi"

source(here::here("scripts/sims-variance-multi/0_utils.R"))

sims <- collect_sim_results(output_name)

sagg <- sims %>%
  group_by(runID, when_run  ) %>%
  summarise(n = n(),
            n_meth = n_distinct(method),
            t_time = sum( time_secs ) / 60,
            .groups = "drop")
sagg

# How long did a run take vs when it was run
plt <- ggplot( data = sagg, aes( when_run, t_time ) ) +
  geom_point()
print( plt )


if ( FALSE ) {


  sagg <- sagg %>%
    mutate( trun = t_time / n )
  ggplot( data = sagg, aes( when_run, trun ) ) +
    geom_point()

  table( sims$overlap_label, sims$tx_type, useNA = "always" )

  n_distinct(sims$runID)

  table( table( sims$runID) )



  # How many rows per run vs when it was run
  ggplot( data = sagg, aes( when_run, n ) ) +
    geom_point()

  filter( sagg, n < 300 )  %>%
    arrange( desc( when_run ) )

  filter( sims, runID == 1 ) %>%
    dplyr::select( runID:method, time_secs )

  qplot( sagg$t_time )

  ggplot( sagg, aes( runID, t_time ) ) +
    geom_point()

  qplot( sims$time_secs  )

  library( tidyverse )
  ggplot( sagg, aes( runID, n ) ) +
    geom_point()


}
