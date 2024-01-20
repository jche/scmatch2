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

table_3 <- res %>%
  mutate(bias_before = (att_est - att_true),
         bias_after = (att_est_debiased - att_true)) %>%
  group_by(name, mu_model) %>%
  summarise(mean_se_boot = mean(sd_boot),
            se_se_boot = ,
            mean_bias_before = mean(bias_before),
            se_bias_before = MCSE_bias(bias_before),
            mean_bias_after = mean(bias_after),
            se_bias_after = MCSE_bias(bias_after),
            true_se = sd(att_est_debiased),
            se_true_se = SE_SEtrue(att_est_debiased) )


## Select residual bootstrap data, and calculate the table 3



# Why is there downward bias in Kang true?
kang_true <- res %>%
  filter(name == "kang", mu_model=="kang_correct")

hist(kang_true$att_est)
