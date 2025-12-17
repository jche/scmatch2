
source( here::here( "scripts/ferman-analysis/02-core-ferman-analysis.R" ) )

library(haven)

options(list(dplyr.summarise.inform = FALSE))
theme_set( theme_classic() )


## Number of used controls
d <- result_table(
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
