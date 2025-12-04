
# Generate figures for supplemental lalonde analysis


# scripts/lalonde-analysis/03-lalonde-figures.R
library(tidyverse)
library(CSM)
library(latex2exp)
library(here)
library(gridExtra)
# library(wesanderson) # Optional, using default colors per previous instruction

options(list(dplyr.summarise.inform = FALSE))
theme_set(theme_classic())

# Source the core analysis
source(here::here("scripts/lalonde-analysis/02-core-lalonde-analysis.R"))

# Ensure figures directory exists
if(!dir.exists(here::here("figures/"))) dir.create(here::here("figures/"), recursive = TRUE)


summary( lalonde_scm )


# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 1. Love Plot (Covariate Balance) - Figure 1 in Paper ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---

# Map X names back to readable names for the plot
covs_names <- c("Black", "Hispanic", "Married", "No Degree",
                "Age", "Education", "Earnings (1974)", "Earnings (1975)",
                "High School Grad", "College Grad")
names(covs_names) <- c( paste0("X", 1:8), "HS", "COLL" )

p_love <- love_plot(lalonde_scm, covs = c( "X5", "X6", "X7", "X8"),
                    covs_names = covs_names[5:8] ) +
  labs(title = "Covariate Balance (Lalonde)")

p_love

ggsave(here::here("figures/lalonde_love.pdf"), plot = p_love, height = 3, width = 6)



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 1.(b) Look at how distances are distributed after matching ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---

dplt <- distance_density_plot( lalonde_scm ) +
  theme(legend.position = "bottom")

dplt

attr( dplt, "table" )

ggsave(here::here("figures/lalonde_distance_density.pdf"), plot = dplt, width = 6, height = 2)



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 2. Comparison of Methods Table ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---

# The paper compares ESS of SCM, Average, and 1-NN. We need to run the
# other two to get these numbers.

ss <- sensitivity_table( lalonde_scm, outcome="Y",
                         include_distances = TRUE )

ss


plt <- ess_plot( lalonde_scm )
plt

plt <- ess_plot( lalonde_scm, feasible_only = TRUE )
plt



# Demo code: Some extra comparisons and explorations

bad_matches( lalonde_scm, 1.5 )

# If we want, we can explore sets by passing lists of IDs
ss <- get_match_sets( lalonde_scm,
                      c( "U00149", "U00065" ),
                      nonzero_weight_only = FALSE )

#ss <- ss %>%
#  dplyr::select( all_of( c( names(DIST_SCALING ), "Z" ) ) )
ss
gen_dm( ss, covs = names( DIST_SCALING ),
        treatment = "Z",
        scaling = DIST_SCALING,
        metric = METRIC )

# Note: We are using maximum distance metric, hence the frequency of
# ties.



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 3. Distances Histogram (Top K) ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---


p_dist_all <- caliper_distance_plot( lalonde_scm, tops = c(1, 3, 6) )
p_dist_all

ggsave(here::here("figures/lalonde_top_k_distances.pdf"), plot = p_dist_all, width = 7, height = 3.5)





# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 4. SATT / Caliper Trade-off Plot (Figure 2) - Cumulative (FSATT -> SATT) ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---

plot_max_cal <- feasible_plot(lalonde_scm, caliper_plot = "both")
p_satt = plot_max_cal$plot_SATT
plot_max_cal = plot_max_cal$plot_max_caliper_size

p_tradeoff_all <- gridExtra::grid.arrange(plot_max_cal, p_satt, ncol=2)


ggsave(here::here("figures/lalonde_fsatt_tradeoff.pdf"), plot = p_tradeoff_all, width = 10, height = 5)



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 5. Sensitivity to caliper plot ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---


# Note: Bit computationally intensive
plt <- caliper_sensitivity_plot( lalonde_scm, lalonde_df,
                         focus = "FATT",
                         R = 11,
                         min_cal = 0, max_cal = 1 )

plt

ggsave( plt, filename = here::here( "figures/lalonde_caliper_sensitivity_plot.pdf"),
        width = 7, height = 5 )

tbl = caliper_sensitivity_table( plt )
tbl

plt_stats <- caliper_sensitivity_plot_stats( plt )
plt_stats
ggsave( plt_stats, filename = here::here( "figures/lalonde_caliper_sensitivity_plot_stats.pdf"),
        width = 7, height = 5 )




