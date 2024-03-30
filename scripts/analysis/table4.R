setwd("~/Dropbox (Harvard University)/Xiang_Luke/scmatch2")
library(tidyr)
library(dplyr)
source("scripts/analysis/MCSE_functions.R")

# load data
I = 100; B=40
BOOT_MTD = "A-E"
# BOOT_MTD = "regression_debiased"
# BOOT_MTD = "regression"
# BOOT_MTD = "cluster"
# BOOT_MTD = "naive"
# FNAME =
#   paste0("data/outputs/sim_toy_results/",BOOT_MTD,"_toy_low_mid_high.csv")
FNAME =
  paste0("data/outputs/bootstrap_toy_rad_new/",BOOT_MTD,"_toy_low_mid_high.csv")
res <- read.csv(file =  FNAME)

table_4 <- res %>%
  mutate(bias_and_var = (att_est - att_true)) %>%
  group_by(deg_overlap) %>%
  summarise(
    sd_att_true = sd(att_true),
    true_se_full = sd(att_est),
    true_se = sd(bias_and_var),
    se_true_se = MCSE_SEtrue(bias_and_var),
    coverage= mean(covered),
    sd_bias_and_var = sd(bias_and_var),
    mean_bias=mean(bias),
    sd_bias = sd(bias-noise),
    mean_noise = mean(noise),
    sd_noise = sd(noise),
    mean_se_boot = mean(sd_boot))
