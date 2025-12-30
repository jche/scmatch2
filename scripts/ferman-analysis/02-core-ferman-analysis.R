# scripts/ferman-analysis/02-core-ferman-analysis.R
# Conduct the core analysis for prepping summary results


library(CSM)
library(tidyverse)

ferman_for_analysis <-
  readRDS(here::here( "scripts/ferman-analysis/data/ferman_for_analysis.rds" ))


match_covs = c("y2007", "y2008", "y2009", "is_sao_paolo")

# Summarize all the match covariates, calculating mean and sd for each
covs_summary <- ferman_for_analysis %>%
  pivot_longer(all_of(match_covs), names_to = "covariate", values_to = "value") %>%
  group_by(covariate) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            sd   = sd(value,   na.rm = TRUE),
            .groups = "drop")
covs_summary

c <- 0.35
covariate_caliper <- c(rep(0.2, 3), 1/1000)
scaling <- 1/covariate_caliper
scaling

ferman_csm <- ferman_for_analysis %>%
  get_cal_matches(
    covs = match_covs,
    treatment = "Z",
    caliper = c,
    metric = "maximum",   # "maximum", "euclidean", "manhattan"
    rad_method = "adaptive",
    scaling = scaling,
    est_method = "scm")



