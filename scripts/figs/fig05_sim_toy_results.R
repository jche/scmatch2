# Plotting script: scripts/figs/fig_sim_toy.R
# full analysis of simulation results
library(tidyverse)
source("scripts/lib/plot_sim.R")

# res_toy <- readr::read_csv("data/outputs/sims-bias_mse-no-filter-unmatched/toy_combined.csv")
# res_toy <- readr::read_csv("data/outputs/sims-bias_mse-21Nov2025/toy_combined.csv") %>%
#   select(-bal2)
# res_toy_26Apr <- readr::read_csv("data/outputs/sims-bias_mse-26Apr2026/toy_combined.csv") %>%
#   select(-bal2)
res_toy <- readr::read_csv("data/outputs/sims-bias_mse/toy_combined.csv",
                           show_col_types = FALSE)
res_L <- pivot_longer(res_toy,
                      diff:last_col(), names_to="method")

res_L %>% group_by( method ) %>%
  summarise( mean_na = mean( is.na( value ) ),
             SE = sd( value, na.rm=TRUE ),
             sd_tau = sd( true_ATT ) ) %>%
  knitr::kable()




if ( FALSE ) {
  # Look at the distribution of estimates to see if we have outliers
  ggplot( res_L, aes( method, value ) ) +
    geom_boxplot() +
    coord_flip()
}

# Windsorize
iqr = IQR( res_L$value, na.rm=TRUE )
iqr
q25 = quantile(res_L$value, 0.25, na.rm=TRUE)
q25
res_L$value[ res_L$value < q25 - 2*iqr ] <- q25 - 2*iqr


org_df <- summarize_bias_rmse( res_L )
org_df <- org_df %>%
  mutate( group = case_when(
    method == "csm_scm" ~ "CSM",
    method == "causal_forest" ~ "CF",
    method == "cem_avg" ~ "CEM",
  str_detect( method, "ipw" ) ~ "IPW",
  TRUE ~ "Other" ) )
org_df


contour_plot( org_df, add_labels = FALSE )
contour_plot(org_df, add_labels = FALSE, focus = c("csm_scm", "cem_avg"))

# debugonce( RMSE_plot )
RMSE_plot( org_df,
           title = "Toy Example",
           xlab = "Value",
           ylab = "Method",
           legend.position = c(0.75, 0.25),
           table_long = TRUE )


ggsave("figures/sim_toy_results.pdf", width = 3.5, height = 3.5)



########################################
# scratch: diagnostic table ----
########################################
res_toy_L <- res_toy %>%
  pivot_longer(
    cols = diff:last_col(), # From diff to the last column
    names_to = "method",
    values_to = "estimate"
  )
res_toy_L

library( simhelpers )
tmp <- res_toy_L %>%
  group_by(method) %>%
  mutate( estimate = estimate - true_ATT,
          true_ATT = 0 ) %>%
  summarise(
    calc_absolute(
      estimates = estimate, true_param = true_ATT,
      criteria = c("bias","stddev", "rmse")
    )
  ) %>%
  arrange(rmse)
tmp


