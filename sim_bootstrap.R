
# check validity of bootstrap

require(tidyverse)
require(mvtnorm)

require(tictoc)

source("R/distance.R")
source("R/sc.R")
source("R/matching.R")
source("R/estimate.R")
source("R/inference.R")
source("R/sim_data.R")
source("R/wrappers.R")
source("R/utils.R")

# parallelize
if (T) {
  library(foreach)
  library(doParallel)
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
}



# function to get bootstrap results ---------------------------------------

res_fun <- function(d, dist_scaling, name, att0, B=1000) {
  
  # get CSM results
  preds_csm <- get_cal_matches(
    df = d,
    metric = "maximum",
    dist_scaling = dist_scaling,
    cal_method = "fixed",
    est_method = "scm",
    return = "all",
    knn = 25)
  
  # bootstrap FSATT
  boot_fsatt <- attr(preds_csm, "scweights") %>% 
    bind_rows() %>% 
    filter(subclass %in% attr(preds_csm, "feasible_subclasses")) %>% 
    agg_co_units() %>% 
    boot_bayesian(B=B)
  
  if (att0) {
    att <- 0
  } else {
    att <- d %>% 
      filter((!Z) | id %in% attr(preds_csm, "feasible_units")) %>% 
      summarize(att = mean(Y1-Y0)) %>% 
      pull(att)
  }
  
  tibble(
    id = i,
    data = name,
    att = att,
    att_hat = get_att_ests(preds_csm),
    sd = sd(boot_fsatt)
  ) %>% 
    mutate(covered = att_hat-1.96*sd <= att | att_hat+1.96*sd >= att,
           length = 1.96*2*sd)
}


# run simulations ---------------------------------------------------------

# repeatedly run for all combinations of pars
res <- foreach(
  i=1:40, 
  .packages = c("tidyverse", "mvtnorm"),
  .combine=rbind) %dopar% {
    source("R/distance.R")
    source("R/sc.R")
    source("R/matching.R")
    source("R/estimate.R")
    source("R/inference.R")
    source("R/sim_data.R")
    source("R/wrappers.R")
    source("R/utils.R")
    
    # for (i in 1:250) {
    
    # toy example
    df <- gen_df_adv(
      nc=500, 
      nt=100, 
      f0_sd = 0.5,
      tx_effect_fun = function(X1, X2) {3*X1+3*X2},
      f0_fun = function(x,y) {
        matrix(c(x,y), ncol=2) %>%
          dmvnorm(mean = c(0.5,0.5),
                  sigma = matrix(c(1,0.8,0.8,1), nrow=2)) * 20   # multiply for more slope!
      })
    num_bins <- 6
    dist_scaling <- df %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) num_bins / (max(x) - min(x))
                         else 1000
                       }))
    
    res_toy <- res_fun(df, dist_scaling, name="toy", att0=F)
    
    # kang
    df <- gen_df_kang(n = 1000)
    nbins <- 5
    dist_scaling <- df %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) nbins / (max(x) - min(x))
                         else 1000
                       }))
    
    res_kang <- res_fun(df, dist_scaling, name="kang", att0=T)
    
    # hain
    df <- gen_df_hain(
      nc=250, 
      nt=50, 
      sigma_e = "n100",   # high overlap condition
      outcome = "nl2",
      sigma_y = 1,
      ATE = 0)
    
    nbins <- 5
    dist_scaling <- df %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) nbins / (max(x) - min(x))
                         else 1000
                       }))
    
    res_hain <- res_fun(df, dist_scaling, name="hain", att0=T)
    
    # acic
    df <- gen_df_acic(
      model.trt="step", 
      root.trt=0.35, 
      overlap.trt="full",
      model.rsp="step", 
      alignment=0.75, 
      te.hetero="high",
      random.seed=i,        # NOTE: set random seed!
      n=1000, 
      p=10
    )
    
    nbins <- 5
    dist_scaling <- df %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) nbins / (max(x) - min(x))
                         else 1000
                       }))
    
    res_acic <- res_fun(df, dist_scaling, name="acic", att0=F)
    
    res <- bind_rows(res_toy, res_kang, res_hain, res_acic)
    
    FNAME <- "bootstrap_results/bootstrap_results.csv"
    if (file.exists(FNAME)) {
      write_csv(res, FNAME, append=T)
    } else {
      write_csv(res, FNAME)
    }
    
    res
  }

stopCluster(cl)



# analyze results ---------------------------------------------------------





