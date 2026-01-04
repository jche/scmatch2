
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


summary( ferman_csm )


# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 1. Love Plot (Covariate Balance)
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---

# Map X names back to readable names for the plot
covs_names <- c( paste0("X", 1:3, " (200", 7:9, " score)"),
                 "X4 (pct Sao Paolo)")

p_love <- love_plot(ferman_csm,
                    covs = c("y2007", "y2008", "y2009"),
                    covs_names = covs_names[1:3] ) +
  labs(title = "Covariate Balance (Ferman)") +
  theme( legend.position = "none" )

p_love

ggsave(here::here("figures/ferman_love.pdf"), plot = p_love, height = 3, width = 6)



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 1.(b) Look at how distances are distributed after matching ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---

dplt <- distance_density_plot( ferman_csm ) +
  theme(legend.position = "bottom")

dplt

attr( dplt, "table" )

ggsave(here::here("figures/ferman_distance_density.pdf"), plot = dplt, width = 6, height = 2)



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 2. Comparison of Methods Table ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---

ss <- sensitivity_table( ferman_csm, outcome="Y",
                         include_distances = TRUE )

ss


s_fin <- ss %>% dplyr::select( Estimate, ATT, SE ) %>%
  mutate( feasible = ifelse( grepl( "FATT", Estimate ), "Feasible", "All" ),
          method = case_when(
            grepl( "1nn", Estimate ) ~ "1-NN",
            grepl( "raw", Estimate ) ~ "Average",
            TRUE ~ "CSM"
          ) ) %>%
  dplyr::select( feasible, method, ATT, SE ) %>%
  pivot_wider( names_from = c( feasible ),
               values_from = c( ATT, SE ), names_vary = "slowest" )

s_fin

# Table of point estimates in the paper
# TODO: Add in the 1-nn SEs --- need to use some other method to get them than our variance estimator?
library( xtable )

print( xtable( s_fin, digits=3 ), include.rownames = FALSE )



plt <- ess_plot( ferman_csm )
plt

plt <- ess_plot( ferman_csm, feasible_only = TRUE )
plt



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 3. Distances Histogram (Top K) ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---


p_dist_all <- caliper_distance_plot( ferman_csm, tops = c(1, 3, 6) )
p_dist_all

ggsave(here::here("figures/ferman_top_k_distances.pdf"), plot = p_dist_all, width = 7, height = 3.5)





# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 4. SATT / Caliper Trade-off Plot - Cumulative (FSATT -> SATT) ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---

plot_max_cal <- feasible_plot(ferman_csm, caliper_plot = "both")
p_satt = plot_max_cal$plot_SATT
plot_max_cal = plot_max_cal$plot_max_caliper_size

p_tradeoff_all <- gridExtra::grid.arrange(plot_max_cal, p_satt, ncol=2)


ggsave(here::here("figures/ferman_fsatt_tradeoff.pdf"), plot = p_tradeoff_all, width = 10, height = 5)



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 5. Sensitivity to caliper plot ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---


# Note: Bit computationally intensive
sens_tbl <- caliper_sensitivity_table( ferman_csm, ferman_for_analysis,
                                       outcome="Y",
                                 R = 11,
                                 min_cal = 0, max_cal = 1 )


plt <- caliper_sensitivity_plot( sens_tbl,
                                 focus = "FATT",
                                 R = 11,
                                 min_cal = 0, max_cal = 1 )

plt


ggsave( plt, filename = here::here( "figures/ferman_caliper_sensitivity_plot.pdf"),
        width = 7, height = 5 )

tbl = caliper_sensitivity_table( plt )
tbl

plt_stats <- caliper_sensitivity_plot_stats( plt )
plt_stats
ggsave( plt_stats, filename = here::here( "figures/ferman_caliper_sensitivity_plot_stats.pdf"),
        width = 7, height = 5 )




