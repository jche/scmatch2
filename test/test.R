source("test/test_bootstrap.R")
source("test/test_sim_data.R")
test_boot_by_resids()
test_boot_bayesian()

boot_otsu_to_test <-
  boot_CSM(dgp_name="otsu", 
           att0=T,
           I=100,
           B=250,
           mu_model="linear",
           boot_mtd="Bayesian",
           n_split=2)
sd(boot_otsu_to_test$att_est_debiased)
mean(boot_otsu_to_test$sd_boot)
mean(boot_otsu_to_test$covered)
mean(boot_otsu_to_test$upper - boot_otsu_to_test$lower)

test_gen_df_otsu()