
# Code timing check


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

# Scaling setup matches the paper's V matrix:
# V = diag {K, K, K, K, 1/3, 1, 1/5000, 1/5000}, where K=1000.
DIST_SCALING <- tibble(
  X1 = 1000, # % Black (Exact Match)
  X2 = 1000, # % Hispanic (Exact Match)
  X3 = 1000, # % Married (Exact Match)
  X4 = 1000, # % No degree (Exact Match)
  X5 = 1/3,  # Average age (Caliper of 3 years)
  X6 = 1,    # Average years of education (Caliper of 1 year)
  X7 = 1/5000, # Average earnings, 1974 (Caliper of $5,000)
  X8 = 1/5000  # Average earnings, 1975 (Caliper of $5,000)
)

# timing check: does saving the distance matrix help at all?

tictoc::tic()

# Run matching
lalonde_scmA <- lalonde_df %>%
  get_cal_matches(
    form = Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8,
    caliper = CALIPER, # 💥 UPDATED to 1
    metric = METRIC,
    rad_method = CAL_METHOD,
    scaling = DIST_SCALING,
    est_method = "scm"
  )

lalonde_scmA
tictoc::toc()

tictoc::tic()

# Run matching
lalonde_scmB <- lalonde_df %>%
  get_cal_matches(
    form = Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8,
    caliper = CALIPER, # 💥 UPDATED to 1
    metric = METRIC,
    rad_method = CAL_METHOD,
    scaling = DIST_SCALING,
    est_method = "scm",
    dm = lalonde_scm$dm_uncapped
  )

lalonde_scmB
tictoc::toc()

expect_equal( lalonde_scmA, lalonde_scmB )


