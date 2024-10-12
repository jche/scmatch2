library(haven)
library(CSM)
library(latex2exp)
library(tidyverse)

options(list(dplyr.summarise.inform = FALSE))
theme_set( theme_classic() )

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

## Number of used controls
d <- full_unit_table(
  ferman_scm,
  feasible_only = TRUE
)
(n_t_SCM <- length(unique(
  (d %>% filter(Z==0, weights!=0))$id )))
(n_t_avg <-length(unique(
  (d %>% filter(Z==0))$id )))


list_distinct_control_1nn <- d %>%
  group_by(Z,subclass) %>%
  filter(Z == 0) %>%
  filter(dist == min(dist)) %>%
  slice(1) %>%
  mutate(weights = 1) %>%
  ungroup() %>%
  distinct(id)

(n_t_1nn <- nrow(list_distinct_control_1nn))
