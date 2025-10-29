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

full_unit_table(ferman_scm, nonzero_weight_only = TRUE) %>%
  rename(isp = is_sao_paolo)

get_ATT_estimate(ferman_scm, "Z", "y2010")

## Number of used controls
summary(ferman_scm)

######
# Histogram ----
######
feasible_subclasses <- feasible_unit_subclass(ferman_scm)
(n_feasible <- length(feasible_subclasses))

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

# Use caliper_table() function which is the proper way to get adacal
ggd_att <-
  result_table(ferman_scm) %>%
  left_join(caliper_table(ferman_scm), by = c("id", "subclass")) %>%
  group_by(subclass) %>%
  summarize(
    adacal = last(adacal),
    tx = y2010[Z == 1] - mean(y2010[Z == 0])
  ) %>%
  arrange(adacal) %>%
  mutate(
    order = 1:n(),
    cum_avg = cumsum(tx) / order
  ) %>%
  filter(subclass %in% feasible_subclasses) %>%
  arrange(order) %>%
  slice(n_feasible:n())

ggd_att

plot_max_caliper_size <-
  ggd_att %>%
  ggplot(aes(x=order, y=adacal)) +
  geom_line(alpha=0.5) +
  geom_point(size=3) +
  theme_classic() +
  labs(y = "Maximum caliper size used",
       x = "Total number of treated units used") +
  expand_limits(x=0, y=0)

feasible_w_adacal <-
  full_unit_table(ferman_scm) %>%
  left_join(caliper_table(ferman_scm), by = c("id", "subclass")) %>%
  group_by(subclass) %>%
  mutate(adacal = ifelse(is.na(adacal), first(na.omit(adacal)), adacal)) %>%
  ungroup() %>%
  arrange(desc(adacal))

feasible_w_adacal$hat_mu_0 <- 0
n_unique_subclass <- length(unique(feasible_w_adacal$subclass))

# Pre-allocate a numeric vector to store the results
se_AEs <- numeric(n_unique_subclass)

# Iterate over the range from the number of unique subclass levels down to n_feasible
for (i in n_unique_subclass:n_feasible) {
  top_subclasses <- feasible_w_adacal %>%
    group_by(subclass) %>%
    summarise(min_adacal = min(adacal, na.rm = TRUE)) %>%
    arrange(min_adacal) %>%
    slice(1:i) %>%
    pull(subclass)

  df_curr <- feasible_w_adacal %>%
    filter(subclass %in% top_subclasses)

  # Get SE
  source(here::here("R/estimate.R"))

  preds_csm_filtered <- df_curr %>%
    filter(Z == 0) %>%
    group_by(subclass) %>%
    filter(n() >= 2) %>%
    ungroup()

  weighted_var <- function(x, wt) {
    n <- length(x)
    wm <- weighted.mean(x, wt)
    sum(wt * (x - wm)^2) * n / (n-1)
  }

  cluster_var_df <-
    preds_csm_filtered %>%
    group_by(subclass) %>%
    summarise(nj = n(),
              w_nj = ess(weights),
              var_cluster = var(y2010), .groups="drop")

  weighted_var_df <- cluster_var_df %>%
    summarise(weighted_var = weighted.mean(var_cluster, w = w_nj), .groups="drop")
  sigma_hat <- sqrt(weighted_var_df$weighted_var)

  Ns <- calc_N_T_N_C(df_curr)
  sd_curr <- sqrt((1/Ns$N_T + 1/Ns$N_C_tilde)) * sigma_hat

  se_AEs[i] <- sd_curr
}

# Output the results
foo <-
  ggd_att %>%
  ggplot(aes(x=order, y=cum_avg)) +
  geom_line(alpha=0.5) +
  geom_point(size=3) +
  theme_classic() +
  labs(y = "Cumulative ATT Estimate",
       x = "Total number of treated units used") +
  expand_limits(color=1)

plot_SATT <- foo +
  geom_hline(yintercept=0, lty="dotted") +
  theme(
    legend.direction="horizontal",
    legend.position.inside = c(0.5, 0.85),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  ) +
  labs(
    y = "Cumulative ATT Estimate",
    x = "Total number of treated units used"
  ) +
  geom_errorbar(
    aes(
      ymin = cum_avg - 1.96 * se_AEs[n_feasible:n_unique_subclass],
      ymax = cum_avg + 1.96 * se_AEs[n_feasible:n_unique_subclass]
    ),
    width = 0.5,
    linewidth = 1
  ) +
  ylim(c(-0.1, 0.2))

# Combine plots using grid.arrange
plot_all <- gridExtra::grid.arrange(
  plot_max_caliper_size,
  plot_SATT,
  ncol=2
)
