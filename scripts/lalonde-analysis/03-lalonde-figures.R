# scripts/lalonde-analysis/03-lalonde-figures.R
library(CSM)
library(latex2exp)
library(tidyverse)
library(here)
library(gridExtra)
# library(wesanderson) # Optional, using default colors per previous instruction

options(list(dplyr.summarise.inform = FALSE))
theme_set(theme_classic())

# Source the core analysis
source(here::here("scripts/lalonde-analysis/02-core-lalonde-analysis.R"))

# Ensure figures directory exists
if(!dir.exists(here::here("figures/lalonde"))) dir.create(here::here("figures/lalonde"), recursive = TRUE)

# -------------------------------------------------------------------------
# 1. Love Plot (Covariate Balance) - Figure 1 in Paper
# -------------------------------------------------------------------------

# Map X names back to readable names for the plot
covs_names <- c("Black", "Hispanic", "Married", "No Degree",
                "Age", "Education", "Earnings (1974)", "Earnings (1975)")
names(covs_names) <- paste0("X", 1:8)

# Using the CSM love_plot logic
# The paper emphasizes the balance of the *feasible* set.
set.seed(2)
p_love <- love_plot(lalonde_scm, covs = paste0("X", 1:8), B = NA) +
  scale_y_discrete(labels = covs_names) +
  labs(title = "Covariate Balance (Lalonde)")

ggsave(here::here("figures/lalonde/lalonde_love.png"), plot = p_love, height = 4, width = 6)


# -------------------------------------------------------------------------
# 2. Comparison Table (Replicating Table 2 in Paper)
# -------------------------------------------------------------------------
# The paper compares SCM, Average, and 1-NN. We need to run the other two.

# A. SCM (Already run)
ess_scm <- result_table(lalonde_scm, feasible_only=TRUE) %>%
  filter(Z==0) %>% summarise(ess = sum(weights)^2/sum(weights^2)) %>% pull(ess)

# B. Average (Simple average of matches in caliper)
lalonde_avg <- lalonde_df %>%
  get_cal_matches(
    caliper = lalonde_params$caliper,
    metric = lalonde_params$metric,
    rad_method = "adaptive", # Keep same radius logic
    scaling = lalonde_params$dist_scaling,
    est_method = "average" # <--- CHANGE
  )
ess_avg <- result_table(lalonde_avg, feasible_only=TRUE) %>%
  filter(Z==0) %>% summarise(ess = sum(weights)^2/sum(weights^2)) %>% pull(ess)

# C. 1-NN (Nearest Neighbor)
# We use rad_method="1nn" and average weighting (since it's 1 unit, avg=1)
lalonde_1nn <- lalonde_df %>%
  get_cal_matches(
    caliper = lalonde_params$caliper,
    metric = lalonde_params$metric,
    rad_method = "1nn",  # <--- CHANGE
    scaling = lalonde_params$dist_scaling,
    est_method = "average"
  )
ess_1nn <- result_table(lalonde_1nn, feasible_only=TRUE) %>%
  filter(Z==0) %>% summarise(ess = sum(weights)^2/sum(weights^2)) %>% pull(ess)

# Print to console to verify against paper (Expect: ~72.6, ~129.9, ~67.0)
print(tibble(Method = c("SCM", "Average", "1-NN"), ESS = c(ess_scm, ess_avg, ess_1nn)))


# -------------------------------------------------------------------------
# 3. Distances Histogram (Top K)
# -------------------------------------------------------------------------

# Extract distance matrix from slot
dist_matrix <- data.frame(t(as.matrix(lalonde_scm$dm_uncapped)))
dm_col_sorted <- apply(dist_matrix, 2, sort)

# Updated plotting function to show the "hard to match" units
# The paper mentions max adaptive caliper is ~2.62.
plot_dm <- function(dist_to_plot){
  tibble(d = as.numeric(as.matrix(dist_to_plot))) %>%
    filter(d < 10) %>% # Filter out exact match dummies (1000), keep adaptive ones (~2.6)
    ggplot(aes(d)) +
    geom_histogram(color="black", binwidth=0.1) +
    geom_vline(xintercept = lalonde_params$caliper, col="red") +
    theme_classic() +
    labs(y=NULL, x = TeX("$d(X_t, X_j)$")) +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
}

# Top 1, 2, 3 distances
p_dist_1 <- plot_dm(dm_col_sorted[1,])
p_dist_2 <- plot_dm(dm_col_sorted[2,])
p_dist_3 <- plot_dm(dm_col_sorted[3,])

p_dist_all <- gridExtra::grid.arrange(p_dist_1, p_dist_2, p_dist_3, ncol=1)
ggsave(here::here("figures/lalonde/lalonde_top_k_distances.png"), plot = p_dist_all, width = 5, height = 8)


# -------------------------------------------------------------------------
# 4. SATT / Caliper Trade-off Plot (Figure 2) - Cumulative (FSATT -> SATT)
# -------------------------------------------------------------------------

# Define 'ess' function (must be here before the loop)
ess <- function( weights ) {
  sum( weights )^2 / sum( weights^2 )
}

# 1. Prepare the full data and Calculate Cumulative Stats GLOBALLY
# We must calculate order and cum_avg on the FULL set first,
# so that the "first" point in our plot represents the FSATT + 1 hard unit.
ggd_att_full <- result_table(lalonde_scm) %>%
  left_join(lalonde_scm$treatment_table, by = c("subclass")) %>%
  group_by(subclass) %>%
  summarize(
    adacal = last(na.omit(adacal)),
    tx = Y[Z==1] - mean(Y[Z==0])
  ) %>%
  ungroup() %>%
  arrange(adacal) %>% # Important: Sort by difficulty
  mutate(
    # This order represents the total number of treated units used
    order = 1:n(),
    # Cumulative ATT Estimate for the set 1:n
    cum_avg = cumsum(tx) / order
  )

# 2. Identify the "Base" set (Feasible units)
# These units are ALWAYS included in the calculation
base_feasible_subclasses <- ggd_att_full %>%
  filter(adacal <= CALIPER) %>%
  pull(subclass)

# 3. Create the plotting data (The "Hard" units)
# We filter for units > CALIPER, but keep the 'cum_avg' calculated from the full set
ggd_att_filtered <- ggd_att_full %>%
  filter(adacal > CALIPER)

# 4. Prepare for SE Loop
n_filtered <- nrow(ggd_att_filtered)
se_estimate_ATT_filtered <- numeric(n_filtered)

# Full unit table for pulling raw data
all_units_for_se <- result_table(lalonde_scm) %>%
  left_join(lalonde_scm$treatment_table, by = c("subclass","id")) %>%
  group_by(subclass) %>%
  mutate(adacal = ifelse(is.na(adacal), first(na.omit(adacal)), adacal)) %>%
  ungroup() %>%
  arrange(adacal)

# 5. Loop: Calculate SEs using (Base Set + Hard Units 1:i)
for (i in 1:n_filtered) {
# i <- 1
  # A. Identify the specific "Hard" subclasses added at this step
  current_hard_subclasses <- ggd_att_filtered$subclass[1:i]

  # B. COMBINE Base Set + Current Hard Set
  # This ensures i=1 includes ALL feasible units plus the 1st hard unit
  combined_subclasses <- c(base_feasible_subclasses, current_hard_subclasses)

  # C. Filter the raw data to this combined set
  df_curr <- all_units_for_se %>%
    filter(subclass %in% combined_subclasses)

  # Check for valid variance groups (at least one control group with n>=2)
  has_valid_variance_group <- df_curr %>%
    filter(Z == 0) %>%
    group_by(subclass) %>%
    summarise(count = n()) %>%
    filter(count >= 2) %>%
    nrow() > 0

  if (!has_valid_variance_group) {
    se_estimate_ATT_filtered[i] <- NA
    next
  }

  # D. TOTAL VARIANCE ESTIMATION
  # We apply the V estimator to the cumulative dataset
  estimate_ATT_result <- get_estimate_ATT(
    matches = df_curr,
    outcome = "Y",
    treatment = "Z",
    variance_method = "pooled"
  )

  if (is.na(estimate_ATT_result$SE)) {
    se_estimate_ATT_filtered[i] <- NA
  } else {
    se_estimate_ATT_filtered[i] <- estimate_ATT_result$SE
  }
}

# 6. Plotting
plot_max_cal <- ggd_att_filtered %>%
  ggplot(aes(x=order, y=adacal)) +
  geom_line(alpha=0.5) +
  geom_point(size=1) +
  theme_classic() +
  labs(y = "Maximum caliper size used", x = "Total number of treated units used")

p_satt <- ggd_att_filtered %>%
  ggplot(aes(x=order, y=cum_avg)) +
  geom_line(alpha=0.5) +
  geom_point(size=1) +
  geom_hline(yintercept=0, lty="dotted") +
  geom_errorbar(
    aes(
      ymin = cum_avg - 1.96 * se_estimate_ATT_filtered,
      ymax = cum_avg + 1.96 * se_estimate_ATT_filtered
    ),
    width = 0.5, alpha = 0.2
  ) +
  labs(y = "Cumulative ATT Estimate", x = "Total number of treated units used") +
  theme_classic()

p_tradeoff_all <- gridExtra::grid.arrange(plot_max_cal, p_satt, ncol=2)
ggsave(here::here("figures/lalonde/lalonde_fsatt_tradeoff.png"), plot = p_tradeoff_all, width = 10, height = 5)

# # Output final point estimate for the full set (SATT)
# full_set <- result_table(lalonde_scm, feasible_only = FALSE)
# print("SATT Estimate (Using all units, for comparison):")
# print(get_att_point_est(full_set))
# -------------------------------------------------------------------------
# 5. Final numbers for text
# -------------------------------------------------------------------------
# FSATT (Feasible only)
feasible <- result_table(lalonde_scm, feasible_only = TRUE)
print("FSATT Estimate (Matches Paper ~$1595):")
print(get_att_point_est(feasible))

# SATT (Full)
full_set <- result_table(lalonde_scm, feasible_only = FALSE)
print("SATT Estimate (Matches Paper ~$1344):")
print(get_att_point_est(full_set))
