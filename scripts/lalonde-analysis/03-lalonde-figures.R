
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
# 1. Love Plot (Covariate Balance) ----
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
# 1. Distribution of distances after matching ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---

dplt <- distance_density_plot( lalonde_scm ) +
  theme(legend.position = "bottom")

dplt

attr( dplt, "table" )

ggsave(here::here("figures/lalonde_distance_density.pdf"),
       plot = dplt, width = 6, height = 2)



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 2. Comparison of Methods Table ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---

# The paper compares ESS and distance of SCM, Average, and 1-NN.

ss <- sensitivity_table( lalonde_scm,
                         include_distances = TRUE )

ss


plt <- ess_plot( lalonde_scm )
plt

plt <- ess_plot( lalonde_scm, feasible_only = TRUE )
plt



# Demo code: Some extra comparisons and explorations to look at
# specific matched sets
if ( FALSE ) {
  bad_matches( lalonde_scm, 1.5 )

  # If we want, we can explore sets by passing lists of IDs
  ss <- result_table( lalonde_scm,
                        id = c( "U00175", "U00162", "U00065" ),
                        nonzero_weight_only = FALSE )

  #ss <- ss %>%
  #  dplyr::select( all_of( c( names(DIST_SCALING ), "Z" ) ) )
  ss
  a <- gen_dm( ss,
               covs = names( DIST_SCALING ),
               treatment = "Z",
               scaling = DIST_SCALING,
               metric = METRIC )


  ss = filter( result_table(lalonde_scm), id %in% txs ) %>%
    dplyr::select( id:COLL )
  ss2 = ss
  ss2$Z = 0
  ss = bind_rows( ss, ss2 )
  a <- gen_dm( ss,
               covs = names( DIST_SCALING ),
               treatment = "Z",
               scaling = DIST_SCALING,
               metric = METRIC )

  # Note: We are using maximum distance metric, hence the frequency of
  # ties.


  # TODO:
  ss <- ss %>%
    dplyr::select( id:COLL )
  ss
  rm = update_matches( lalonde_scm, data=ss )
  rm
  result_table( rm )
  gen_dm( ss,
          covs = names( DIST_SCALING ),
          treatment = "Z",
          scaling = DIST_SCALING,
          metric = METRIC )
}



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 3. Distances Histogram (Top K) ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---


p_dist_all <- caliper_distance_plot( lalonde_scm, tops = c(1, 3, 6) )
p_dist_all

ggsave(here::here("figures/lalonde_top_k_distances.pdf"),
       plot = p_dist_all, width = 7, height = 3.5)





# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 4. FATT to ATT plot
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---

sensitivity_table( lalonde_scm, outcome="Y" )

plot_max_cal <- feasible_plot(lalonde_scm, caliper_plot = "both")

p_satt = plot_max_cal$plot_SATT
plot_max_cal = plot_max_cal$plot_max_caliper_size

plot_max_cal
p_satt


p_tradeoff_all <- gridExtra::grid.arrange(plot_max_cal, p_satt, ncol=2)
p_tradeoff_all

ggsave(here::here("figures/lalonde_fsatt_tradeoff.pdf"),
       plot = p_tradeoff_all, width = 10, height = 5)



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---
# 5. Sensitivity to caliper plot ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*---


# Note: Bit computationally intensive
cal_sens_tbl <- caliper_sensitivity_table( lalonde_scm, lalonde_df,
                                           R = 30,
                                           outcome = "Y",
                                           min_cal = 0, max_cal = 10.01 )

cal_sens_tbl

plt_cal_sens <- caliper_sensitivity_plot( cal_sens_tbl,
                                          focus = c( "ATT", "FATT" ) ) +
  facet_wrap( ~ feasible ) +
  geom_vline( xintercept = max( lalonde_scm$adacalipers ), lty=3, col="black" )


plt_cal_sens


# Alt plot
caliper_sensitivity_plot( cal_sens_tbl ) +
  geom_vline( xintercept = max( lalonde_scm$adacalipers ), lty=3, col="black" )



plt_cal_sens



a = get_distance_table( lalonde_scm )
mean( a$closest <= 1 )

max( a$closest )

ggsave( plt_cal_sens, filename = here::here( "figures/lalonde_caliper_sensitivity_plot.pdf"),
        width = 7, height = 5 )


plt_stats <- caliper_sensitivity_plot_stats( cal_sens_tbl )
plt_stats
ggsave( plt_stats, filename = here::here( "figures/lalonde_caliper_sensitivity_plot_stats.pdf"),
        width = 7, height = 5 )




# Addendum: Looking at control unit reuse ----


ctls <- result_table( lalonde_scm, return = "agg_co_units" ) %>%
  filter( Z == 0 )

summary( ctls$weights )

ggplot( ctls, aes( weights ) ) +
  geom_histogram( binwidth = 0.5 ) +
  scale_x_continuous( breaks = 0:10 )
nrow(ctls)
sum( ctls$weights >= 1 )

sum( ctls$weights )^2 / sum( ctls$weights^2 )

estimate_ATT( lalonde_scm )

rr = result_table( lalonde_scm )
table( table( rr$id ) )

summary( lalonde_scm )


# Addendum: Looking at bad matches and balance for those matches ----


lalonde_scm

bad_tx <- bad_matches( lalonde_scm, threshold = 1 )
bad_tx
sum( bad_tx$Z )

bad_tx %>%
  pivot_longer( cols = all_of( names(DIST_SCALING) ),
                names_to="Variable",
                values_to="Value" ) %>%
  group_by( Variable, Z ) %>%
  summarise( Mean = weighted.mean( Value, weights ),
             SD = sqrt( weighted.mean( (Value - Mean)^2,
                                       weights ) ),
             N = sum( weights ),
             .groups = "drop" ) %>%
  arrange( Variable, Z ) %>%
  pivot_wider( names_from = Z,
               values_from = c( Mean, SD, N ),
               names_prefix = "Z_" ) %>%
  mutate( Diff_in_Means = Mean_Z_1 - Mean_Z_0 )


# Looking at impact as a function of caliper ----

plt_imp <- impact_curve( lalonde_scm, outcome="Y" )
plt_imp

imp_tbl <- attr( plt_imp, "table" )
imp_tbl

summary( imp_tbl$precision )
