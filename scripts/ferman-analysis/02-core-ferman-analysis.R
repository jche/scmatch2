
# Conduct the core analysis for prepping summary results


library(CSM)
library(tidyverse)

ferman_for_analysis <-
  readRDS(here::here( "scripts/ferman-analysis/data/ferman_for_analysis.rds" ))

c <- 0.35
covariate_caliper <- c(rep(0.2, 3), 1/1000)
scaling <- 1/covariate_caliper

ferman_scm <- ferman_for_analysis %>%
  get_cal_matches(
    covs = c("y2007", "y2008", "y2009", "is_sao_paolo"),
    treatment = "Z",
    caliper = c,
    metric = "maximum",   # "maximum", "euclidean", "manhattan"
    rad_method = "adaptive",
    scaling = 1/covariate_caliper,
    est_method = "scm")
