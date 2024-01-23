setwd("~/Dropbox (Harvard University)/Xiang_Luke/scmatch2")
library(tidyverse)

# load data
I = 40; B=40
FNAME =
  paste0("./sim_toy_results/kang_toy_naive_I_",I,"_B_",B,".csv")
res <- read.csv(file =  FNAME)


kurtosis <- function(x){
  S_T = sd(x)
  kurt = mean( (x - mean(x))^4 ) / S_T^4
  return(kurt)
}

MCvar_SEtrue <- function(x){
  S_T = sd(x); R <- length(x); k_T <- kurtosis(x)
  return(S_T^2 * sqrt( (k_T-1)/R ))
}

MCSE_SEtrue <- function(x){
  return(sqrt(MCvar_SEtrue(x)))
}

MCSE_bias <- function(x){
  S_T = sd(x); R <- length(x)
  return( sqrt(S_T^2 /R ))
}

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
# kang_correct <- kangs %>% filter(mu_model=="kang_correct")
# hist(kang_correct$att_est_debiased)
# kang_wrong <- kangs %>% filter(mu_model=="linear")
# hist(kang_wrong$att_est_debiased)

toys <- res_res_boot %>%
  filter(name=="toy", n_split==1, boot_mtd=="wild")

# # Assess variation of true means of toys
# hist(toys$att_true)
# mean(toys$att_true)
# sd(toys$att_true)
# hist(toys$att_est)
# sd(toys$att_est)
#
# mean(toys$att_est - 3)
# sd((toys$att_est - 3))
# mean(toys$att_est - toys$att_true)
# sd((toys$att_est - toys$att_true))


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
