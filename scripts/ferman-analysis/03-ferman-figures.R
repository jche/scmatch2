
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
if(!dir.exists(here::here("figures/"))) {
  dir.create(here::here("figures/"), recursive = TRUE)
}

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
# 1. Distance distribution after matching ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

dplt <- distance_density_plot( ferman_csm ) +
  theme(legend.position = "bottom")

dplt

attr( dplt, "table" )

ggsave(here::here("figures/ferman_distance_density.pdf"),
       plot = dplt, width = 6, height = 2)



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# 2. Comparison of Methods Table ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

summary( ferman_csm, outcome="Y" )

estimate_ATT( ferman_csm )

ferman2 <- update_matches( ferman_csm, data=ferman_for_analysis, k = 2 )
ferman2


ss <- sensitivity_table( ferman_csm, outcome="Y",
                         include_distances = TRUE )

ss

ss2 <- sensitivity_table( ferman2, outcome="Y",
                          include_distances = TRUE )
ss2$Estimate = paste0( ss2$Estimate, " (k=2)" )

# ESS table
ss %>%
  bind_rows(  ss2 ) %>%
  dplyr::select( Estimate, mean_dist, median_dist, N_T, N_C, ESS_C, p_drop ) %>%
  arrange( Estimate ) %>%
  knitr::kable( digits = 3 )



ferman_csm

s_fin <- ss %>%
  bind_rows( ss2[c(1,4),] ) %>%
  dplyr::select( Estimate, ATT, SE ) %>%
  mutate( feasible = ifelse( grepl( "FATT", Estimate ),
                             "Feasible", "All" ),
          method = case_when(
            grepl( "1nn", Estimate ) ~ "1-NN",
            grepl( "raw", Estimate ) ~ "Average",
            grepl( "k=2", Estimate ) ~ "CSM (k=2)",
            TRUE ~ "CSM"
          ) ) %>%
  dplyr::select( feasible, method, ATT, SE ) %>%
  pivot_wider( names_from = c( feasible ),
               values_from = c( ATT, SE ), names_vary = "slowest" )

# Initial table of point estimates
s_fin


# Now add SEs for 1-nn

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

# compare to our version
ferman_1nn <- get_cal_matches( ferman_for_analysis,
    covs = match_covs,
    treatment = "Z",
    caliper = c,
    metric = "euclidean",
    rad_method = "knn",
    scaling = scaling,
    k = 1,
    est_method = "average")
estimate_ATT( ferman_1nn )
s_fin$ATT_All[ s_fin$method == "1-NN" ]


matched_pairs <- tibble(
  treated = ferman_for_analysis$id[m_out$index.treated],
  control = ferman_for_analysis$id[ m_out$index.control]
)
matched_pairs
fm = result_table( ferman_1nn ) %>%
  dplyr::select( id, Z, subclass ) %>%
  pivot_wider( names_from = Z, values_from=id )
fm

fmm = left_join( fm, matched_pairs, by = c("1" = "treated") )
mean( fmm$`0` == fmm$control )
dels <- filter( fmm, `0` != control )
idd = c( dels$`1`, dels$`0`, dels$control )
filter( ferman_for_analysis, id %in% as.character( dels[1,] ) )
filter( ferman_for_analysis, id %in% as.character( dels[2,] ) )


# Now subset to feasible units and then do nearest neighbor
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


# compare to our version
ferman_1nn <- get_cal_matches( ferman_for_analysis,
                               covs = match_covs,
                               treatment = "Z",
                               caliper = c,
                               metric = "maximum",
                               rad_method = "knn-capped",
                               scaling = scaling,
                               k = 1,
                               est_method = "scm")
estimate_ATT( ferman_1nn )
s_fin$ATT_All[ s_fin$method == "1-NN" ]



m_out2$est
m_out2$se.standard


# Plug in the SEs
s_fin$ATT_All[ s_fin$method == "1-NN" ] = m_out$est
s_fin$SE_All[ s_fin$method == "1-NN" ] = m_out$se.standard
s_fin$ATT_All[ s_fin$method == "1-NN" ] - m_out2$est
s_fin$ATT_Feasible[ s_fin$method == "1-NN" ] = m_out2$est
s_fin$SE_Feasible[ s_fin$method == "1-NN" ] = m_out2$se.standard


## Table of treatment impacts by method ----

s_fin

# Make latex form
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

ggsave(here::here("figures/ferman_top_k_distances.pdf"),
       plot = p_dist_all, width = 7, height = 3.5)





# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# 4. SATT / Caliper Trade-off Plot - Cumulative (FSATT -> SATT) ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

plot_max_cal <- feasible_plot(ferman_csm, caliper_plot = "both")
p_satt = plot_max_cal$plot_SATT
plot_max_cal = plot_max_cal$plot_max_caliper_size

p_tradeoff_all <- gridExtra::grid.arrange(plot_max_cal, p_satt, ncol=2)


ggsave(here::here("figures/ferman_fsatt_tradeoff.pdf"),
       plot = p_tradeoff_all, width = 10, height = 5)



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
ggsave( plt_stats,
        filename = here::here( "figures/ferman_caliper_sensitivity_plot_stats.pdf"),
        width = 7, height = 5 )





# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# 6. Comparing to CEM ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-


source( here::here( "scripts/lib/wrappers.R" ) )


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
cem


ferman_cem <- ferman_for_analysis %>%
  get_cal_matches(
    covs = match_covs,
    treatment = "Z",
    caliper = c,
    metric = "maximum",   # "maximum", "euclidean", "manhattan"
    rad_method = "adaptive",
    scaling = scaling,
    est_method = "scm")



# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Histogram of number of controls used by treated units ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

feasible_subclasses <- feasible_unit_subclass(ferman_csm)
n_feasible <- length(feasible_subclasses)
n_feasible

matched_controls <- result_table(ferman_csm, "all") %>%
  group_by(subclass) %>%
  summarise(n_controls = n()-1) %>%
  filter(subclass %in% feasible_subclasses)

sum(matched_controls$n_controls)
hist_matched_controls <-
  ggplot(matched_controls, aes(x = n_controls)) +
  geom_histogram(binwidth = 3, fill = "skyblue", color = "white") +
  labs(
    title = "Number of controls used by treated units",
    x = "Number of Matched Controls",
    y = "Frequency"
  ) +
  theme_minimal()

matched_controls_nonzero_weights <-
  ferman_csm$matches %>%
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
    title = "Number post synthetic control step",
    x = "Number of Matched Controls",
    y = "Frequency"
  ) +
  theme_minimal()

p <- gridExtra::grid.arrange(
  hist_matched_controls,
  hist_matched_controls_nonzero_weights,
  ncol=2
)
p

ggsave(
  here::here("figures/hist-n-co.pdf"),
  plot=p,
  width=8.1,
  height=5.3
)





# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Table: Number of used controls ----
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-


library(haven)

options(list(dplyr.summarise.inform = FALSE))
theme_set( theme_classic() )


## Number of used controls
d <- result_table(
  ferman_csm,
  feasible_only = TRUE
)
d

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

