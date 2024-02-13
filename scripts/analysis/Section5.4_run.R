setwd("~/Dropbox (Harvard University)/Xiang_Luke/scmatch2")
library(tidyverse)
require(mvtnorm)
source("R/distance.R")
source("R/sc.R")
source("R/matching.R")
source("R/estimate.R")
source("R/inference.R")
source("R/sim_data.R")
source("R/wrappers.R")
source("R/utils.R")
source("R/bootstrap.R")

save_res_to_csv<-
  function(curr_res,
           FNAME){
    if (file.exists(FNAME)) {
      write_csv(curr_res, FNAME, append=T)
    } else {
      write_csv(curr_res, FNAME)
    }
  } # save_res_to_csv

I = 100; B=40;
BOOT_MTD = "A-E"
# BOOT_MTD = "regression_debiased"
# BOOT_MTD = "regression"
# BOOT_MTD = "cluster"
# BOOT_MTD = "naive"
FNAME =
  paste0("./sim_toy_results/",BOOT_MTD,"_toy_low_mid_high.csv")
source("R/bootstrap.R")
source("R/sim_data.R")
toy_naive_low<-
  boot_CSM(dgp_name="toy",
           att0=F,
           I=I,
           B=B,
           mu_model="linear",
           boot_mtd=BOOT_MTD,
           n_split=1,
           kang_true = F,
           toy_ctr_dist = 0.5)
save_res_to_csv(toy_naive_low,
                FNAME = FNAME)

toy_naive_mid<-
  boot_CSM(dgp_name="toy",
           att0=F,
           I=I,
           B=B,
           mu_model="linear",
           boot_mtd=BOOT_MTD,
           n_split=1,
           kang_true = F,
           toy_ctr_dist = 0.3)
save_res_to_csv(toy_naive_mid,
                FNAME = FNAME)

toy_naive_high<-
  boot_CSM(dgp_name="toy",
           att0=F,
           I=I,
           B=B,
           mu_model="linear",
           boot_mtd=BOOT_MTD,
           n_split=1,
           kang_true = F,
           toy_ctr_dist = 0.1)
save_res_to_csv(toy_naive_high,
                FNAME = FNAME)



I = 40; B=40
FNAME =
  paste0("./sim_toy_results/kang_toy_naive_I_",I,"_B_",B,".csv")
source("R/bootstrap.R")
toy_naive<-
  boot_CSM(dgp_name="toy",
           att0=F,
           I=1,
           B=2,
           mu_model="linear",
           boot_mtd="cluster",
           n_split=1,
           kang_true = F)
save_res_to_csv(toy_naive,
                FNAME = FNAME)


source("R/bootstrap.R")
kang_naive<-
  boot_CSM(dgp_name="kang",
           att0=T,
           I=I,
           B=B,
           mu_model="kang_correct",
           boot_mtd="naive",
           n_split=1,
           kang_true = T)
save_res_to_csv(
  kang_naive,
  FNAME =FNAME)


source("R/bootstrap.R")
kang_false_naive<-
  boot_CSM(dgp_name="kang",
           att0=T,
           I=I,
           B=B,
           mu_model="linear",
           boot_mtd="naive",
           n_split=1,
           kang_true = F)
save_res_to_csv(
  kang_false_naive,
  FNAME = FNAME)


# Get results using
I <- 100; B=20
FNAME <- paste0("./sim_toy_results/toy_bootstrap.csv")
source("R/bootstrap.R")
kang_false_resid<-
  boot_CSM(dgp_name="kang",
           att0=T,
           I=I,
           B=B,
           mu_model="linear",
           boot_mtd="wild",
           n_split=1,
           kang_true = F)
save_res_to_csv(
  kang_false_resid,
  FNAME = FNAME)

