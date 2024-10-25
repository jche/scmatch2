
#setwd("~/Dropbox (Harvard University)/Xiang_Luke/scmatch2")
library(tidyr)
library(dplyr)
source( here::here( "scripts/analysis/MCSE_functions.R") )

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
R = 100
# FNAME =
#   here::here(
#     paste0("data/outputs/bootstrap_toy_rad_new/",
#            "test_",
#            "A_E_toy_low_mid_high_R=",R,".csv")
#   )
R = 500
FNAME =
  here::here(
    paste0("data/outputs/A-E-overlap-by-prop-unif/",
           "A_E_toy_low_mid_high_R=",R,".csv")
  )
res <- read.csv(file =  FNAME)

table_4 <- res %>%
  mutate(
    bias_and_var = (att_est - att_true),
    deg_overlap = factor(deg_overlap, levels = c("low", "mid", "high"))
  ) %>%
  group_by(deg_overlap) %>%
  summarise(
    SE_est = mean(se_AE),
    SE_True = true_SE[1],
    N_C_tilde = N_C_tilde[1],
    mean_bias = mean(bias),
    coverage = mean(covered),
    coverage_with_true_SE = mean(covered_with_true_SE)
  ) %>%
  arrange(deg_overlap)

table_4

# let plot
res_low <- res %>% filter(deg_overlap=="low")
hist(res_low$att_est)
abline(v=res_low$att_true[1],col="red")
abline(v=mean(res_low$att_est),col="blue")

mean_est = mean(res_low$att_est)
mean_bias = res_low$att_true[1] -mean_est
true_se = sd(res_low$error)
se_est = mean(res_low$se_AE)

pnorm(mean_est + 1.96 * true_se + mean_bias,
      mean=mean_est, sd=true_se ) -
  pnorm(mean_est - 1.96 * true_se + mean_bias,
        mean=mean_est, sd=true_se)
