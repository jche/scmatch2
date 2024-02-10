setwd("~/Dropbox (Harvard University)/Xiang_Luke/scmatch2")
library(tidyr)
library(dplyr)
source("analysis/MCSE_functions.R")

# load data
I = 100; B=40
BOOT_MTD = "A-E"
# BOOT_MTD = "regression_debiased"
# BOOT_MTD = "regression"
# BOOT_MTD = "cluster"
# BOOT_MTD = "naive"
FNAME =
  paste0("./sim_toy_results/",BOOT_MTD,"_toy_low_mid_high.csv")
res <- read.csv(file =  FNAME)
res$deg_overlap <- c(rep("low", 100),
                     rep("medium", 100),
                     rep("high", 100))

table_4 <- res %>%
  mutate(bias_before = (att_est - att_true),
         bias_after = (att_est_debiased - att_true)) %>%
  group_by(deg_overlap) %>%
  summarise(mean_se_boot = mean(sd_boot),
            true_se = sd(att_est_debiased),
            se_true_se = MCSE_SEtrue(att_est_debiased))

