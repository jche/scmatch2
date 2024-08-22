


library(tidyverse)
require(mvtnorm)

library( CSM )
source( "scripts/analysis/boot_CSM_simulation_code.R" )

save_res_to_csv<-
  function(curr_res,
           FNAME){
    if (file.exists(FNAME)) {
      write_csv(curr_res, FNAME, append=TRUE)
    } else {
      write_csv(curr_res, FNAME)
    }
  } # save_res_to_csv


I = 1000; B=40;
BOOT_MTD = "A-E"
# BOOT_MTD = "regression_debiased"
# BOOT_MTD = "regression"
# BOOT_MTD = "cluster"
# BOOT_MTD = "naive"
FNAME =
  here::here( paste0("data/outputs/bootstrap_toy_rad_new/",BOOT_MTD,"_toy_low_mid_high.csv") )


sample_dat <- get_df_scaling_from_dgp_name(dgp_name="toy",
                             kang_true=FALSE,
                             toy_ctr_dist=0.5)
names(sample_dat)
samp_dat <- sample_dat$df_dgp
samp_dat
ggplot( samp_dat, aes( X1, X2, col=as.factor(Z) ) ) +
  geom_point() +
  coord_fixed()
samp_dat$tau = samp_dat$Y1 - samp_dat$Y0
skimr::skim( samp_dat )

# Test simulation driver ----

if ( FALSE ) {
  tst <-
    boot_CSM(dgp_name="toy",
             att0=F,
             I=5,
             B=10,
             mu_model="linear",
             boot_mtd="A-E",
             n_split=1,
             kang_true = FALSE,
             toy_ctr_dist = 0.5)

}

# Run the simulation for inference ----

run_full_simulation <- function( dgp_name ) {

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
toy_naive_low$deg_overlap <- "low"
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
toy_naive_mid$deg_overlap <- "mid"
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
toy_naive_high$deg_overlap <- "high"
save_res_to_csv(toy_naive_high,
                FNAME = FNAME)


}


run_full_simulation( "toy" )
run_full_simulation( "kang" )


####
# Inference on Kang-Schafer ----

#   The results of this DGP are no longer in the main paper anymore.
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

kang_false_resid<-
  boot_CSM(dgp_name="kang",
           att0=T,
           I=I,
           B=B,
           mu_model="linear",
           boot_mtd="wild",
           n_split=1,
           kang_true = F)

FNAME <- here::here( paste0("./sim_toy_results/toy_bootstrap.csv") )
save_res_to_csv(
  kang_false_resid,
  FNAME = FNAME)

