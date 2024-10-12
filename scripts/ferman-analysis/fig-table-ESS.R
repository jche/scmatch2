library(CSM)
library(latex2exp)
library(tidyverse)

ferman_for_analysis <-
  readRDS(here::here( "data/inputs/ferman_for_analysis.rds" ))

c <- 0.35
covariate_caliper <- c(rep(0.2, 3), 1/1000)
ferman_scm <- ferman_for_analysis %>%
  get_cal_matches(
    covs = c("y2007", "y2008", "y2009", "is_sao_paolo"),
    treatment = "Z",
    caliper = c,
    metric = "maximum",   # "maximum", "euclidean", "manhattan"
    rad_method = "adaptive",
    scaling = 1/covariate_caliper,
    est_method = "scm",
    return = "sc_units")

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
