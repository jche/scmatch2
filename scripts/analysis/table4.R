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
# FNAME =
#   paste0("data/outputs/bootstrap_toy_rad_new/",BOOT_MTD,"_toy_low_mid_high.csv")
FNAME =
  here::here(
    paste0("data/outputs/bootstrap_toy_rad_new/",
           "test_",
           "A_E_toy_low_mid_high.csv")
  )
res <- read.csv(file =  FNAME)

table_4 <- res %>%
  mutate(bias_and_var = (att_est - att_true)) %>%
  group_by(deg_overlap) %>%
  summarise(
    RMSE = sd(bias_and_var),
    SE_est = mean(se_AE),
    sd_noise = sd(error),
    mean_bias=mean(bias),
    coverage= mean(covered),
    # mean_bias_sq = mean(bias^2),
    # sd_att_true = sd(att_true),
    # true_se_full = sd(att_est),
    # se_true_se = MCSE_SEtrue(bias_and_var),

    # sd_bias_and_var = sd(bias_and_var),
    # sd_bias = sd(bias-noise),
    # mean_noise = mean(noise)
    )

# let plot
res_low <- res %>% filter(deg_overlap=="low")
hist(res_low$att_est)
abline(v=res_low$att_true[1],col="red")
abline(v=mean(res_low$att_est),col="blue")

# i think the next is to simulate this 1000 times
# mean while i play, i do the BG stuff. i do weekly review.
# blank tomorrow afternoon out for coding learning
