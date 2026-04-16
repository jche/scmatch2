
library(haven)
library(CSM)
library(latex2exp)
library(tidyverse)

options(list(dplyr.summarise.inform = FALSE))
theme_set(theme_classic())

# Source the core analysis instead of duplicating code
source(here::here("scripts/ferman-analysis/02-core-ferman-analysis.R"))

# Drop unneeded variables
ferman_for_analysis$UF = NULL
ferman_for_analysis$y2011 = NULL
ferman_for_analysis$y2012 = NULL
ferman_for_analysis$Control = NULL
ferman_for_analysis$Y = NULL

names(ferman_for_analysis)

# Display main results
ferman_scm

result_table(ferman_scm, nonzero_weight_only = TRUE) %>%
  rename(isp = is_sao_paolo)

estimate_ATT(ferman_scm, treatment="Z", outcome="y2010")

## Number of used controls
summary(ferman_scm)


######
# Histogram ----
######

feasible_subclasses <- feasible_unit_subclass(ferman_scm)
n_feasible <- length(feasible_subclasses)
n_feasible

matched_controls <- result_table(ferman_scm, "all") %>%
  group_by(subclass) %>%
  summarise(n_controls = n()-1) %>%
  filter(subclass %in% feasible_subclasses)

sum(matched_controls$n_controls)
hist_matched_controls <-
  ggplot(matched_controls, aes(x = n_controls)) +
  geom_histogram(binwidth = 3, fill = "skyblue", color = "white") +
  labs(
    x = "Number of Matched Controls",
    y = "Frequency"
  ) +
  theme_minimal()

matched_controls_nonzero_weights <-
  ferman_scm$matches %>%
  bind_rows() %>%
  filter(weights > 0) %>%
  group_by(subclass) %>%
  summarise(n_controls = n()-1) %>%
  filter(subclass %in% feasible_subclasses)
sum(matched_controls_nonzero_weights$n_controls)

hist_matched_controls_nonzero_weights <-
  ggplot(matched_controls_nonzero_weights, aes(x = n_controls)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "white") +
  labs(
    x = "Number of Matched Controls",
    y = "Frequency"
  ) +
  theme_minimal()

p <- gridExtra::grid.arrange(
  hist_matched_controls,
  hist_matched_controls_nonzero_weights,
  ncol=2
)
ggsave(
  here::here("figures/hist-n-co.png"),
  plot=p,
  width=8.1,
  height=5.3
)


#####
# SATT plot ----
######

plot_max_caliper_size <-
  feasible_plot( ferman_scm, caliper_plot=TRUE )

plot_SATT <-
  feasible_plot( ferman_scm, caliper_plot=FALSE )

# Combine plots using grid.arrange
plot_all <- gridExtra::grid.arrange(
  plot_max_caliper_size,
  plot_SATT,
  ncol=2
)



