
# Generate figures for supplemental ferman analysis


# scripts/ferman-analysis/03-ferman-figures.R
library(tidyverse)
library(CSM)
library(latex2exp)
library(here)
library(gridExtra)
# library(wesanderson) # Optional, using default colors per previous instruction

options(list(dplyr.summarise.inform = FALSE))
theme_set(theme_classic())

# Source the core analysis
source(here::here("scripts/ferman-analysis/02-core-ferman-analysis.R"))

# Ensure figures directory exists
if(!dir.exists(here::here("figures/"))) dir.create(here::here("figures/"), recursive = TRUE)


summary( ferman_scm )


# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 1. Love Plot (Covariate Balance)
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---

# Map X names back to readable names for the plot
covs_names <- c( paste0("X", 1:3, " (200", 7:9, " score)"),
                 "X4 (pct Sao Paolo)")

p_love <- love_plot(ferman_scm,
                    covs = c("y2007", "y2008", "y2009", "is_sao_paolo"),
                    covs_names = covs_names ) +
  labs(title = "Covariate Balance (Ferman)")

p_love

ggsave(here::here("figures/ferman_love.pdf"), plot = p_love, height = 3, width = 6)



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 1.(b) Look at how distances are distributed after matching ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---

dplt <- distance_density_plot( ferman_scm ) +
  theme(legend.position = "bottom")

dplt

attr( dplt, "table" )

ggsave(here::here("figures/ferman_distance_density.pdf"), plot = dplt, width = 6, height = 2)



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 2. Comparison of Methods Table ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---

ss <- sensitivity_table( ferman_scm, outcome="Y",
                         include_distances = TRUE )

ss


plt <- ess_plot( ferman_scm )
plt

plt <- ess_plot( ferman_scm, feasible_only = TRUE )
plt



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 3. Distances Histogram (Top K) ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---


p_dist_all <- caliper_distance_plot( ferman_scm, tops = c(1, 3, 6) )
p_dist_all

ggsave(here::here("figures/ferman_top_k_distances.pdf"), plot = p_dist_all, width = 7, height = 3.5)





# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 4. SATT / Caliper Trade-off Plot - Cumulative (FSATT -> SATT) ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---

plot_max_cal <- feasible_plot(ferman_scm, caliper_plot = "both")
p_satt = plot_max_cal$plot_SATT
plot_max_cal = plot_max_cal$plot_max_caliper_size

p_tradeoff_all <- gridExtra::grid.arrange(plot_max_cal, p_satt, ncol=2)


ggsave(here::here("figures/ferman_fsatt_tradeoff.pdf"), plot = p_tradeoff_all, width = 10, height = 5)



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 5. Sensitivity to caliper plot ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---


# Note: Bit computationally intensive
plt <- caliper_sensitivity_plot( ferman_scm, ferman_for_analysis,
                                 focus = "FATT",
                                 R = 11,
                                 min_cal = 0, max_cal = 1 )

plt


# bug: caliper_sensitivity_table() --> sensitivity_table() calls
#  estimate_ATT(feasible_only=T), which fails if nothing is feasible



ggsave( plt, filename = here::here( "figures/ferman_caliper_sensitivity_plot.pdf"),
        width = 7, height = 5 )

tbl = caliper_sensitivity_table( plt )
tbl

plt_stats <- caliper_sensitivity_plot_stats( plt )
plt_stats
ggsave( plt_stats, filename = here::here( "figures/ferman_caliper_sensitivity_plot_stats.pdf"),
        width = 7, height = 5 )




