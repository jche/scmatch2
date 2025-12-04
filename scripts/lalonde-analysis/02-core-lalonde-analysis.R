# scripts/lalonde-analysis/02-core-lalonde-analysis.R
# Conduct the core analysis for prepping summary results

library(CSM)
library(tidyverse)
library(here)

lalonde_df <- readRDS(here::here("scripts/lalonde-analysis/data/lalonde_for_analysis.rds"))
dim( lalonde_df )

# MATCHING SETTINGS BASED ON LALONDE PAPER SNIPPET (Section 3.2)

# TODO: If we set caliper to 0.5, then we actually see a change in the
CALIPER <- 0.25
METRIC <- "maximum"
CAL_METHOD <- "adaptive" # The paper uses adaptive caliper: c_t = max{1, d_t}

lalonde_df$HS = ifelse( lalonde_df$X6 >= 12, 1, 0 )
lalonde_df$COLL = ifelse( lalonde_df$X6 >= 16, 1, 0 )

# Scaling setup matches the paper's V matrix:
# V = diag {K, K, K, K, 1/3, 1, 1/5000, 1/5000}, where K=1000.
DIST_SCALING <- tibble(
  X1 = 1000, # % Black (Exact Match)
  X2 = 1000, # % Hispanic (Exact Match)
  X3 = 1000, # % Married (Exact Match)
  X4 = 1000, # % No degree (Exact Match)
  X5 = 1/3,  # Average age (Caliper of 3 years)
  X6 = 1/2,    # Average years of education (Caliper of 1 year)
  X7 = 1/5000, # Average earnings, 1974 (Caliper of $5,000)
  X8 = 1/5000,  # Average earnings, 1975 (Caliper of $5,000)
  HS = 1/2,
  COLL = 1/2,
)


# Run matching
lalonde_scm <- lalonde_df %>%
  get_cal_matches(
    form = Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + HS + COLL,
    caliper = CALIPER,
    metric = METRIC,
    rad_method = CAL_METHOD,
    scaling = DIST_SCALING,
    est_method = "scm"
  )

lalonde_scm

cdp <- caliper_distance_plot( lalonde_scm, target_percentile = 0.67 )
cdp

# based on plot
CALIPER = 1/2

lalonde_scm <- update_matches( lalonde_scm, lalonde_df, caliper = CALIPER )
lalonde_scm

