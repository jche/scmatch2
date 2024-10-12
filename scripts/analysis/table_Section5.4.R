

library(tidyr)
library(dplyr)
source( here::here( "scripts/analysis/MCSE_functions.R") )

# load data
R = 1000
FNAME =
  here::here(
    paste0("data/outputs/A-E-overlap-by-prop-unif/",
           "A_E_toy_low_mid_high_R=",R,".csv")
  )
res <- read.csv(file = FNAME)
head( res )

table_section_5_4 <- res %>%
  mutate(
    bias_and_var = (att_est - att_true),
    deg_overlap = factor(deg_overlap, levels = c("low", "mid", "high"))
  ) %>%
  group_by(deg_overlap) %>%
  summarise(
    SE_est = mean(se_AE),
    SE_True = true_SE[1], # Note: N_C_tilde is fixed for all runs, so is true_SE
    N_C_tilde = N_C_tilde[1],
    mean_bias = mean(bias),
    coverage = mean(covered),
    coverage_with_true_SE = mean(covered_with_true_SE)
  ) %>%
  arrange(deg_overlap)

table_section_5_4
