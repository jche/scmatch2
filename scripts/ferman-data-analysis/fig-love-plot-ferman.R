library(latex2exp)
library(tidyverse)
library(CSM)

ferman_for_analysis <-
  readRDS(here::here( "data/inputs/ferman_for_analysis.rds" ))

c <- 0.35
covariate_caliper <- c(rep(0.2, 3), 1/1000)
scaling <- 1/covariate_caliper

ferman_for_analysis <- ferman_for_analysis %>%
  mutate(`X1 (2007 score)` = y2007,
         `X2 (2008 score)` = y2008,
         `X3 (2009 score)` = y2009,
         `X4 (pct Sao Paolo)` = is_sao_paolo)
covs <-c("X1 (2007 score)", "X2 (2008 score)",
         "X3 (2009 score)", "X4 (pct Sao Paolo)")

ferman_scm <- ferman_for_analysis %>%
  get_cal_matches(
    covs = covs,
    treatment = "Z",
    caliper = c,
    metric = "maximum",   # "maximum", "euclidean", "manhattan"
    rad_method = "adaptive",
    scaling = scaling,
    est_method = "scm",
    return = "sc_units")

ferman_scm$id <-
  gsub("*_syn","",ferman_scm$id)
ferman_scm$id <- as.integer(ferman_scm$id)

love_plot(
  res = ferman_scm,
  covs = covs,
  B=NA) +
  ylim(c(-0.0065, 0.006))

ggsave(
  filename="figures/love-plot-ferman.png",
  width=8.76,
  height=5.33
)
