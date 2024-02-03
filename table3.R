setwd("~/Dropbox (Harvard University)/Xiang_Luke/scmatch2")
library(tidyverse)
source("analysis/MCSE_functions.R")

# load data
I = 40; B=40
FNAME =
  paste0("./sim_toy_results/kang_toy_naive_I_",I,"_B_",B,".csv")
res <- read.csv(file =  FNAME)


table_3_naive <- res %>%
  mutate(bias_before = (att_est - att_true),
         bias_after = (att_est_debiased - att_true)) %>%
  group_by(name, mu_model) %>%
  summarise(mean_se_boot_naive = mean(sd_boot))


## Select residual bootstrap data, i.e., toy, kang, kang_true in
#   and calculate the table 3
res_res_boot <-
  read.csv(file =  paste0("./sim_toy_results/toy_bootstrap.csv"))

kangs <- res_res_boot %>%
  filter(name == "kang") %>%
  filter(boot_mtd == "Bayesian" |
           boot_mtd == "wild") %>%
  filter(mu_model=="kang_correct"|
           mu_model=="linear")

toys <- res_res_boot %>%
  filter(name=="toy", n_split==1, boot_mtd=="wild")


table_3_residual <- rbind(kangs,toys) %>%
  mutate(bias_before = (att_est - att_true),
         bias_after = (att_est_debiased - att_true)) %>%
  group_by(name, mu_model) %>%
  summarise(mean_se_boot_residual = mean(sd_boot),
            true_se = sd(att_est_debiased),
            se_true_se = MCSE_SEtrue(att_est_debiased),
            true_se_undebiased = sd(att_est),
            mean_bias_before = mean(bias_before),
            se_bias_before = MCSE_bias(bias_before),
            mean_bias_after = mean(bias_after),
            se_bias_after = MCSE_bias(bias_after),)
table_3 <-
  left_join(table_3_naive,
            table_3_residual,
            by=c("name", "mu_model"))
table_3 <- table_3 %>%
  mutate(across(where(is.numeric), ~signif(., digits = 3))) %>%
  arrange(desc(name),desc(mu_model)) %>%
  mutate(RB =mean_se_boot_naive / true_se)
