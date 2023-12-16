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




# res_bayesian_boot_toy<-boot_bayesian(dgp_name="toy", att0=F,I=100,B=250)
# res_bayesian_boot_kang<-
#   boot_bayesian(dgp_name="kang", att0=T,I=100,B=250)
# res_bayesian_boot_kang_nonlinear<-
#   boot_bayesian(dgp_name="kang", att0=T,I=100,B=250,mu_model="nonlinear")

res_bayesian_boot_kang_correct<-
  boot_bayesian(dgp_name="kang", att0=T,I=100,B=250,mu_model="kang_correct")
res_bayesian_boot_kang_correct$boot_mtd <-
  "Bayesian_mu_correct"
print(mean(res_bayesian_boot_kang_correct$covered))

FNAME = "./sim_toy_results/toy_bootstrap_with_boot_sd.csv"
if (file.exists(FNAME)) {
  write_csv(res_bayesian_boot_kang_correct, FNAME, append=T)
} else {
  write_csv(res_bayesian_boot_kang_correct, FNAME)
}







# note: standard, bayesian, and wild bootstraps all capture true sd
#  - but wild bootstrap requires vec-mean(vec) 
#    BEFORE multiplying by the weight
#  - others work with vec-mean(vec), but it's not required



## Next: make the naive boostrap fail for toy example
I <- 10 # do 100 times
B <- 100
covered <- CI_lower <- CI_upper <- 
  att_true <- att_est <- att_debiased<-numeric(I)
T_boots <- numeric(B)
set.seed(123)
t0 <- proc.time()
for (i in 1:I){
  # i <- 1
  print(i)
  # Generate data of size 100 according to the DGP N(mu,1)
  df_toy <- gen_df_adv(
    nc=500, 
    nt=100, 
    f0_sd = 0.5,
    tx_effect_fun = function(X1, X2) {3*X1+3*X2},
    f0_fun = function(x,y) {
      matrix(c(x,y), ncol=2) %>%
        dmvnorm(mean = c(0.5,0.5),
                sigma = matrix(c(1,0.8,0.8,1), nrow=2)) * 20   # multiply for more slope!
    })
  dist_scaling_toy <- df_toy %>%
    summarize(across(starts_with("X"),
                     function(x) {
                       if (is.numeric(x)) 6 / (max(x) - min(x))
                       else 1000
                     }))
  
  
  # Get the matches 
  preds_csm <- get_cal_matches(
    df = df_toy,
    metric = "maximum",
    dist_scaling = dist_scaling_toy,
    cal_method = "fixed",
    est_method = "scm",
    # return = "agg_co_units",
    return = "all",
    knn = 25)   # NOTE: knn does nothing when cal_method = "fixed"...
  
  att_est[i] <- get_att_ests(preds_csm)
  
  
  # Perform naive bootstrap
  ## Naive bootstrap: Draw 600 from all vs draw 100 trt, 500 control
  for (b in 1:B){
    # b <- 1
    print(b)
    ## Drawing 100 trt, 500 control from df_toy
    # Splitting the dataframe into two subsets
    df_toy_T <- df_toy[df_toy$Z == T, ]
    df_toy_F <- df_toy[df_toy$Z == F, ]
    
    # Sampling with replacement
    sampled_T <- df_toy_T[sample(nrow(df_toy_T), 100, replace = TRUE), ]
    sampled_F <- df_toy_F[sample(nrow(df_toy_F), 500, replace = TRUE), ]
    
    # Combining the samples to create the bootstrapped dataset
    df_boot <- rbind(sampled_T, sampled_F)
    
    # Re-match, and compute the att
    dist_scaling_toy_boot <- df_boot %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) 6 / (max(x) - min(x))
                         else 1000
                       }))
    
    preds_csm <- get_cal_matches(
      df = df_boot,
      metric = "maximum",
      dist_scaling = dist_scaling_toy_boot,
      cal_method = "fixed",
      est_method = "scm",
      # return = "agg_co_units",
      return = "all",
      knn = 25)   
    
    T_boots[b] <- get_att_ests(preds_csm)
    
  }
  
  CI_lower[i] = quantile(T_boots, 0.025)
  CI_upper[i] = quantile(T_boots, 0.975)
  
  # Get true ATT for coverage
  att_true[i] <- att <- df_toy %>%
    # filter((!Z) | id %in% attr(preds_csm, "feasible_units")) %>%
    filter(Z & !(id %in% attr(preds_csm, "unmatched_units"))) %>%
    summarize(att = mean(Y1-Y0)) %>%
    pull(att)
  
  print(paste0("ATT is ", signif(att,3),
               "; LB is ", signif(CI_lower[i],3),
               "; UB is ", signif(CI_upper[i],3)))
  covered[i] = (CI_lower[i] < att) & (att < CI_upper[i])
  
  print(paste0("Covered is ", covered[i]))
  
}
t1 <- proc.time()

res_save_naive_boot <- tibble(id=1:I,
                                 name="toy",
                                 boot_mtd = "Naive",
                                 att_true = att_true,
                                 att_est= att_est,
                                 att_est_debiased = att_debiased,
                                 lower=CI_lower,upper=CI_upper, 
                                 covered=covered)


FNAME = "./sim_toy_results/toy_bootstrap.csv"
if (file.exists(FNAME)) {
  write_csv(res_save_naive_boot, FNAME, append=T)
} else {
  write_csv(res_save_naive_boot, FNAME)
}





# DGP simple data
# input: N, mu
# output: 

# boot simple data 
# input: 
# output: a vector of bootstrapped residuals

N <- 100; B <- 1000; mu <- 1

I <- 100 # do 100 times
covered <- CI_lower <- CI_upper <- numeric(I)
T_star <- numeric(B)
colnames_res <- c("id","data","true_param",
                  "est","q025","q975","covered")
p <- length(colnames_res)
res_local <- array(dim = c(I, p))
set.seed(123)
# input: 
for (i in 1:I){
  print(i)
  # Generate data of size 100 according to the DGP N(mu,1)
  X <- rnorm(N, mean=mu) 
  
  # Get the estimate 
  X_bar = mean(X)
  
  # Get the residuals
  X_resid = X - X_bar
  # Perform Bayesian bootstrap
  for (b in 1:B){
    W = gtools::rdirichlet(1, alpha=rep(1,N))
    T_star[b] = sum(X_resid * W)
  }
  
  CI_lower[i] = X_bar - quantile(T_star, 0.975)
  CI_upper[i] = X_bar - quantile(T_star, 0.025)
  covered[i] = (CI_lower[i] < mu) & (mu < CI_upper[i])
  # return tibble
  # save to CSV
  # if save to local is true, then append to the local df
}
mean(covered)
idx_not_covered<-which(covered!=T)
CI_lower[idx_not_covered]
CI_upper[idx_not_covered]




# goal: estimate mean
vec <- rnorm(N, mean=mu)

# true sd:
map_dbl(
  1:B,
  function(x) {
    rnorm(N) %>% mean(mean=mu)
  }
) %>% 
  sd()

# naive bootstrap
map_dbl(
  1:B,
  function(x) {
    sample(vec, N, replace=T) %>% mean()
  }) %>% 
  sd()

# non-par bootstrap, which is equivalent to naive bootstrap
map_dbl(
  1:B,
  function(x) {
    (vec * as.numeric(rmultinom(1, size=N, prob=rep(1/N,N))) / N) %>% sum()
  }) %>% 
  sd()

# bayesian bootstrap (Rubin, 1981)
map_dbl(
  1:B,
  function(x) {
    (vec * as.numeric(gtools::rdirichlet(1, alpha=rep(1,N)))) %>% sum()
  }) %>% 
  sd()

# wild bootstrap
#  - note: needs to subtract mean to be accurate!
map_dbl(
  1:B,
  function(x) {
    ((vec-mean(vec)) * sample(
      c( -(sqrt(5)-1)/2, (sqrt(5)+1)/2 ),
      prob = c( (sqrt(5)+1)/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5)) ),
      replace = T, size = N)) %>% mean()
  }) %>% 
  sd()

