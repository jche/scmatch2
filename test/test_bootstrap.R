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


test_split_data <- function(){
  dgp_obj <-
    get_df_scaling_from_dgp_name(dgp_name="toy")
  list2env(dgp_obj, envir = environment())


  # split_data
  #   input: df_to_split, n_split
  #   output: a df, df_to_split with one more column of
  #       group label of one of 1,2,...,n_split


  preds_csm <- get_cal_matches(
    df = df_dgp,
    metric = "maximum",
    dist_scaling = dist_scaling,
    cal_method = "fixed",
    est_method = "scm",
    return = "all",
    knn = 25)





  df_dgp_nested <- df_dgp_splitted %>%
    group_by(group_label) %>%
    nest()
  get_SL_fit_nested <- function(){
    # input: df_dgp_nested,
    #         X_names,
    #         Y_name,
    #         SL.library
    # output:
    n.df <- length(dat)
    res <- vector("list", length = n.df)
    for (i in 1:n.df) {
      df <- dat[[i]]
      res[[i]] <- get_SL_fit(df_to_fit,
                             X_names,
                             Y_name,
                             SL.library)
    }
  }
  length(df_dgp_nested)
  df_dgp_nested$mu_fit <-
    pi_fitter(dat_nested$data)

}


test_regression_se <- function(){

  # Get true se of the non-debiased estimators
  #   in kang, which is around 1.81
  FNAME = "./sim_toy_results/toy_bootstrap.csv"
  res_boot_all <-
    read_csv(FNAME)
  res_bayesian_boot_kang_correct <-
    res_boot_all %>% filter(name=="kang",
                            boot_mtd == "Bayesian_mu_correct")
  (sd_att_est<-
      sd(res_bayesian_boot_kang_correct$att_est))
  (sd_att_est_debiased<-
      sd(res_bayesian_boot_kang_correct$att_est_debiased))
  # Generate a kang df
  dgp_obj <-
    get_df_scaling_from_dgp_name(dgp_name="kang")
  list2env(dgp_obj, envir = environment())

  # plot the debiased amount in kang
  hist(res_bayesian_boot_kang_correct$att_est -
         res_bayesian_boot_kang_correct$att_est_debiased)

  #
  library(estimatr)
  lm_test<-
    lm_robust(
      Y ~ Z + (V1 + V2 + V3 + V4),
      # Y ~ Z + (X1 + X2 + X3 + X4),
              df_dgp,
              se_type = "HC0") # Default variance estimator is HC2 robust standard errors
  summary(lm_test)
  preds_csm <- get_cal_matches(
    df = df_dgp,
    metric = "maximum",
    dist_scaling = dist_scaling,
    cal_method = "fixed",
    est_method = "scm",
    return = "agg_co_units",
    knn = 25)
  lm_test_weighted<-
    lm_robust(
      Y ~ Z * (V1 + V2 + V3 + V4),
      # Y ~ Z + (X1 + X2 + X3 + X4),
              data=preds_csm,
              weights=weights,
              se_type = "HC2")
  summary(lm_test_weighted)
  sd_att_est_debiased
  sd_att_est

  # Test on toy
  res_bayesian_boot_toy_correct <-
    res_boot_all %>% filter(name=="toy",
                            boot_mtd=="Bayesian")
  (sd_att_est<-
      sd(res_bayesian_boot_toy_correct$att_est))
  (sd_att_est_debiased<-
    sd(res_bayesian_boot_toy_correct$att_est_debiased))

  # plot the debiased amount in toy
  hist(res_bayesian_boot_toy_correct$att_est -
         res_bayesian_boot_toy_correct$att_est_debiased)
  # Generate a toy df
  dgp_obj <-
    get_df_scaling_from_dgp_name(dgp_name="toy")
  list2env(dgp_obj, envir = environment())

  library(estimatr)
  preds_csm <- get_cal_matches(
    df = df_dgp,
    metric = "maximum",
    dist_scaling = dist_scaling,
    cal_method = "fixed",
    est_method = "scm",
    return = "agg_co_units",
    knn = 25)
  lm_test<-
    lm_robust(Y ~ Z * (X1 + X2 ),
              data=preds_csm,
              # weights=weights,
              se_type = "HC2")
  summary(lm_test)
  lm_test_weighted<-
    lm_robust(Y ~ Z * (X1 + X2 ),
              data=preds_csm,
              weights=weights,
              se_type = "HC2")
  summary(lm_test_weighted)
}




test_SL_fit_and_pred_linear <- function(){
  dgp_obj <-
    get_df_scaling_from_dgp_name(dgp_name="toy")
  list2env(dgp_obj, envir = environment())

  preds_csm <- get_cal_matches(
    df = df_dgp,
    metric = "maximum",
    dist_scaling = dist_scaling,
    cal_method = "fixed",
    est_method = "scm",
    return = "all",
    knn = 25)

  df_controls <- df_dgp[df_dgp$Z == 0,]
  X_names_toy<- c("X1","X2")



  SL_fit_lm <- get_SL_fit(df_to_fit=df_controls,
                       X_names = X_names_toy,
                       Y_name = "Y",
                      SL.library = "SL.lm")



  SL_pred <- get_SL_pred(SL_fit=SL_fit_lm,
                         df_pred=preds_csm,
                         X_names=X_names_toy)
  # Test predictions using lm
  # Using controls to train, and predict for both
  #   control and trts
  lm_0 <- lm(Y ~ X1 + X2,
             data = df_controls)
  test_preds <-
    predict(lm_0,
            newdata = preds_csm[,X_names_toy])
  # Test: the predictions should be the same
  stopifnot(sum(SL_pred[,,drop=T] - test_preds)==0)
}

test_SL_fit_and_pred_non_linear <- function(){
  dgp_obj <-
    get_df_scaling_from_dgp_name(dgp_name="toy")
  list2env(dgp_obj, envir = environment())

  preds_csm <- get_cal_matches(
    df = df_dgp,
    metric = "maximum",
    dist_scaling = dist_scaling,
    cal_method = "fixed",
    est_method = "scm",
    return = "all",
    knn = 25)

  df_controls <- df_dgp[df_dgp$Z == 0,]
  X_names_toy<- c("X1","X2")
  # SL_lib <- c("SL.glm",
  #             "SL.glmnet",
  #             "SL.randomForest",
  #             "SL.xgboost")


  SL_fit <- get_SL_fit(df_to_fit=df_controls,
                          X_names = X_names_toy,
                          Y_name = "Y",
                          SL.library = SL_lib)


  SL_pred <- get_SL_pred(SL_fit=SL_fit,
                         df_pred=preds_csm,
                         X_names=X_names_toy)
  # Test predictions using lm
  # # Using controls to train, and predict for both
  # #   control and trts
  # lm_0 <- lm(Y ~ X1 + X2,
  #            data = df_controls)
  # test_preds <-
  #   predict(lm_0,
  #           newdata = preds_csm[,X_names_toy])
  # # Test: the predictions should be the same
  # stopifnot(sum(SL_pred[,,drop=T] - test_preds)==0)
}


test_boot_wild <- function(){
  wild_boot_toy_to_test <-
    boot_wild(dgp_name="toy",
              att0=F,
              I=100,
              B=200)
  # Test: coverage should be around 95%
  mean(wild_boot_toy_to_test$covered)
}


# source("R/bootstrap.R")
# test_boot_bayesian()

test_boot_bayesian <- function(){
  bayesian_boot_toy_to_test <-
    boot_bayesian(dgp_name="toy",
                  att0=F,
                  I=1,
                  B=250)
  # write_csv(bayesian_boot_toy_to_test,
  #           file="./test/test_data/boot_toy_test.csv")
  bayesian_boot_toy_expected <-
    read_csv(file="./test/test_data/boot_toy_test.csv")

  bayesian_boot_kang_to_test <-
    boot_bayesian(dgp_name="kang",
                  att0=T,
                  I=1,
                  B=250,
                  mu_model = "kang_correct")
  # write_csv(bayesian_boot_kang_to_test,
  #           file="./test/test_data/boot_kang_test.csv")
  bayesian_boot_kang_expected <-
    read_csv(file="./test/test_data/boot_kang_test.csv")

  # Test: both the ATT estimate, and
  #     the upper bound on CI on the first row
  #     should equal
  stopifnot(bayesian_boot_toy_to_test$att_est[1] ==
              bayesian_boot_toy_expected$att_est[1])
  stopifnot(bayesian_boot_toy_to_test$upper[1] ==
              bayesian_boot_toy_expected$upper[1])
  stopifnot(bayesian_boot_kang_to_test$att_est[1] ==
              bayesian_boot_kang_expected$att_est[1])
  stopifnot(bayesian_boot_kang_to_test$upper[1] ==
              bayesian_boot_kang_expected$upper[1])
}

test_boot_by_resids<-function(){
  # Test 1: Gaussian
  sd_X = 2; n = 500
  I <- 100
  mean_X <- se_boot <- numeric(I)
  for (i in 1:I){
    X <- rnorm(n, mean = 0, sd= sd_X)
    mean_X[i] <- mean(X)
    T_star <- boot_by_resids(resids=X-mean(X),
                             B=250,
                             boot_mtd="Bayesian",
                             seed_addition=i)
    se_boot[i] = sd(T_star)
  }
  sd(mean_X)

  se_gaussian_expected <- sd / sqrt(n)
  # Test: se should be around sd / sqrt(n)
  tol <- 0.01
  stopifnot(
    abs(mean(se_boot) - se_gaussian_expected)
            <
              tol)
  print("Test boot by resids passed")
}

