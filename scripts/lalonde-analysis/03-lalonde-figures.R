
# Generate figures for supplemental lalonde analysis


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

p_love

ggsave(here::here("figures/lalonde/lalonde_love.png"), plot = p_love, height = 4, width = 6)


# -------------------------------------------------------------------------
# 2. Comparison Table (Replicating Table 2 in Paper)
# -------------------------------------------------------------------------
# The paper compares SCM, Average, and 1-NN. We need to run the other two.

summary( lalonde_scm )
estimate_ATT( lalonde_scm )

tbl <- compare_ess( lalonde_scm, lalonde_df )
print( tbl )




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
dists <- dm_col_sorted[1:3,] %>%
  t() %>%
  as_tibble( .name_repair = "unique") %>%
  set_names( c("Top_1", "Top_2", "Top_3") ) %>%
  pivot_longer( cols=everything(),
                names_to = "Rank",
                values_to = "Distance")

summary(dists$Distance)

ggplot( dists, aes( Distance )  ) +
  facet_wrap( ~Rank, ncol=1 ) +
  geom_histogram( color="black" )



# Old style histogram
# TODO/NOTE: I don't like the shifting x-axis.  Replace with the above
p_dist_1 <- plot_dm(dm_col_sorted[1,])
p_dist_2 <- plot_dm(dm_col_sorted[2,])
p_dist_3 <- plot_dm(dm_col_sorted[3,])

p_dist_all <- gridExtra::grid.arrange(p_dist_1, p_dist_2, p_dist_3, ncol=1)
ggsave(here::here("figures/lalonde/lalonde_top_k_distances.png"), plot = p_dist_all, width = 5, height = 8)


# -------------------------------------------------------------------------
# 4. SATT / Caliper Trade-off Plot (Figure 2) - Cumulative (FSATT -> SATT)
# -------------------------------------------------------------------------

plot_max_cal <- feasible_plot(lalonde_scm, caliper_plot = "both")
p_satt = plot_max_cal$plot_SATT
plot_max_cal = plot_max_cal$plot_max_caliper_size

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
