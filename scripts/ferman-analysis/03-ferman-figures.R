
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


# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# 1. Love Plot (Covariate Balance) ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

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



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# 1.(b) Look at how distances are distributed after matching ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

dplt <- distance_density_plot( ferman_csm ) +
  theme(legend.position = "bottom")

dplt

attr( dplt, "table" )

ggsave(here::here("figures/ferman_distance_density.pdf"), plot = dplt, width = 6, height = 2)



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# 2. Comparison of Methods Table ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

ss <- sensitivity_table( ferman_csm, outcome="Y",
                         include_distances = TRUE )

ss

ferman_csm

s_fin <- ss %>% dplyr::select( Estimate, ATT, SE ) %>%
  mutate( feasible = ifelse( grepl( "FATT", Estimate ),
                             "Feasible", "All" ),
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

# First calculate SEs for 1-nn
library(Matching)

X <- as.matrix(ferman_for_analysis[, match_covs])
head( X )
apply( X, 2, sd )
X <- scale(X)

calipers = covariate_caliper * 0.35
calipers

weight.matrix = diag( c( 1, 1, 1, 1000 ) )
weight.matrix

m_out <- Matching::Match(
  Y = ferman_for_analysis$Y,
  Tr = ferman_for_analysis$Z,
  X = X,
  M = 1,
  Weight = 3,
  Weight.matrix =  weight.matrix,
  replace = FALSE,
  ties = FALSE,
  estimand = "ATT"
)
summary( m_out )

# Now do feasible
drop = bad_matches( ferman_csm, 0.35 ) %>%
  pull( subclass ) %>%
  unique()
drop
drop = ferman_for_analysis$id %in% drop
table( drop )

m_out2 <- Matching::Match(
  Y = ferman_for_analysis$Y[!drop],
  Tr = ferman_for_analysis$Z[!drop],
  X = X[!drop,],
  M = 1,
  Weight = 3,
  Weight.matrix =  weight.matrix,
  replace = FALSE,
  ties = FALSE,
  estimand = "ATT"
)
summary( m_out2 )

m_out$est
m_out$se.standard

s_fin$ATT_All[ s_fin$method == "1-NN" ] = m_out$est
s_fin$SE_All[ s_fin$method == "1-NN" ] = m_out$se.standard
s_fin$ATT_Feasible[ s_fin$method == "1-NN" ] = m_out2$est
s_fin$SE_Feasible[ s_fin$method == "1-NN" ] = m_out2$se.standard


library( xtable )

print( xtable( s_fin, digits=3 ), include.rownames = FALSE )


plt <- ess_plot( ferman_csm )
plt

plt <- ess_plot( ferman_csm, feasible_only = TRUE )
plt



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# 3. Distances Histogram (Top K) ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-


p_dist_all <- caliper_distance_plot( ferman_csm, tops = c(1, 3, 6) )
p_dist_all

ggsave(here::here("figures/ferman_top_k_distances.pdf"), plot = p_dist_all, width = 7, height = 3.5)





# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# 4. SATT / Caliper Trade-off Plot - Cumulative (FSATT -> SATT) ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

plot_max_cal <- feasible_plot(ferman_csm, caliper_plot = "both")
p_satt = plot_max_cal$plot_SATT
plot_max_cal = plot_max_cal$plot_max_caliper_size

p_tradeoff_all <- gridExtra::grid.arrange(plot_max_cal, p_satt, ncol=2)


ggsave(here::here("figures/ferman_fsatt_tradeoff.pdf"), plot = p_tradeoff_all, width = 10, height = 5)



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# 5. Sensitivity to caliper plot ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-


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





# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# 6. Comparing to CEM ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-


source( here::here( "scripts/wrappers.R" ) )


scaling
match_covs
cut2007 = seq( min( ferman_for_analysis$y2007 ) - 0.1,
               max( ferman_for_analysis$y2007 )+ 0.1, by = 0.35 * 5 )
cut2008 = seq( min( ferman_for_analysis$y2008 )- 0.1,
               max( ferman_for_analysis$y2008 )+ 0.1, by = 0.35 * 5 )
cut2009 = seq( min( ferman_for_analysis$y2009 )- 0.1,
               max( ferman_for_analysis$y2009 )+ 0.1, by = 0.35 * 5 )
cutISP = c(-0.5,0.5,1.5)
cutpoints = list( y2007 = cut2007,
                  y2008 = cut2008,
                  y2009 = cut2009,
                  is_sao_paolo = cutISP )
cutpoints

cem <- get_cem_matches( data = ferman_for_analysis,
                        covs = match_covs,
                        cutpoints = cutpoints,
                        Z_FORMULA = "Z",
                        est_method = "average",
                        return = "all" )


ferman_cem <- ferman_for_analysis %>%
  get_cal_matches(
    covs = match_covs,
    treatment = "Z",
    caliper = c,
    metric = "maximum",   # "maximum", "euclidean", "manhattan"
    rad_method = "adaptive",
    scaling = scaling,
    est_method = "scm")




