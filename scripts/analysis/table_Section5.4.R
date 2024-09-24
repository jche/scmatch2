setwd("~/Dropbox (Harvard University)/Xiang_Luke/scmatch2")
library(tidyr)
library(dplyr)
source("scripts/analysis/MCSE_functions.R")

# load data
R = 500
FNAME =
  here::here(
    paste0("data/outputs/A-E-overlap-by-prop-unif/",
           "A_E_toy_low_mid_high_R=",R,".csv")
  )
res <- read.csv(file =  FNAME)

table_section_5_4 <- res %>%
  mutate(
    bias_and_var = (att_est - att_true),
    deg_overlap = factor(
      deg_overlap,
      # levels = c("low", "mid", "high")
      levels =  c("very_low", "low", "mid", "high", "very_high")
      )
  ) %>%
  group_by(deg_overlap) %>%
  summarise(
    SE_est = mean(se_AE),
    MCSE_SE_est = sd(se_AE)/n(),
    SE_True = mean(true_SE),
    MCSE_SE_True = sd(true_SE)/n(),
    # MCSE = sd(att_est),
    N_C_tilde = N_C_tilde[1],
    mean_bias = mean(bias),
    coverage = mean(covered),
    coverage_with_true_SE = mean(covered_with_true_SE)
  ) %>%
  arrange(deg_overlap)

write.csv(table_section_5_4,
          file = "tables/table_section_5_4.csv")
