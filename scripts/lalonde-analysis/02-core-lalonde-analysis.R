# scripts/lalonde-analysis/02-core-lalonde-analysis.R
# Conduct the core analysis for prepping summary results

library(CSM)
library(tidyverse)
library(here)

lalonde_df <- readRDS(here::here("scripts/lalonde-analysis/data/lalonde_for_analysis.rds"))
dim( lalonde_df )

# MATCHING SETTINGS BASED ON LALONDE PAPER SNIPPET (Section 3.2)
# The paper uses a fixed initial caliper of c = 1.
CALIPER <- 1
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

# Run matching
lalonde_scm <- lalonde_df %>%
  get_cal_matches(
    caliper = CALIPER, # 💥 UPDATED to 1
    metric = METRIC,
    rad_method = CAL_METHOD,
    scaling = DIST_SCALING,
    est_method = "scm"
  )

lalonde_scm


# Helper to expose the scaling and metric to downstream scripts if needed
lalonde_params <- params( lalonde_scm )

lalonde_params


cc <- caliper_table( lalonde_scm ) %>%
  arrange(-adacal)
cc

bad_matches( lalonde_scm, 10000 )

bad_matches( lalonde_scm, 10000, nonzero_weight_only = FALSE ) %>%
  group_by( subclass ) %>%
  summarise (n = n() )

get_matches( lalonde_scm, "U00149", nonzero_weight_only = FALSE )

result_table( lalonde_scm, include_caliper = TRUE )
%>%
  arrange( desc( caliper ) )
