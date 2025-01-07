library(latex2exp)

source( here::here( "scripts/ferman-analysis/02-core-ferman-analysis.R" ) )

# ESS plots ----
scweights_df <- full_unit_table( ferman_scm )

subclass_feasible <-
  feasible_unit_subclass( ferman_scm )

feasible <-
  full_unit_table(
    ferman_scm,
    feasible_only = TRUE
  )

ess_plot(feasible)
