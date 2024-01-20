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

I = 40; B=40
FNAME = 
  paste0("./sim_toy_results/kang_toy_naive_I_",I,"_B_",B,".csv")
source("R/bootstrap.R")
toy_naive<-
  boot_CSM(dgp_name="toy",
           att0=F,
           I=I,
           B=B,
           mu_model="linear",
           boot_mtd="naive",
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

