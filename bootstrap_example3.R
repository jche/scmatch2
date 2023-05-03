
# back to using the "SC units" version of the bootstrap
#  - now these elements have mean 0 again, so this should work


require(tidyverse)
require(mvtnorm)

source("R/distance.R")
source("R/sc.R")
source("R/matching.R")
source("R/estimate.R")
source("R/inference.R")
source("R/sim_data.R")
source("R/wrappers.R")
source("R/utils.R")

# true sd: 0.123
# wild bootstrap estimated sd: 0.161 (slightly conservative, but ok!)
# regular bootstrap estimated sd: 0.170

# NOTE: when f0_sd is pumped up, estimator becomes anticonservative!
#  - e.g., f0_sd=5: wild bootstrap gives estimate of 0.60, true of 1.34
#     - note: other weighted bootstrap gives estimate of 1.1ish!
dgp <- function() {
  gen_df_adv(
    nc=500,
    nt=100,
    f0_sd = 0.5,
    tx_effect_fun = function(X1, X2) {3*X1+3*X2},
    f0_fun = function(x,y) {
      matrix(c(x,y), ncol=2) %>%
        dmvnorm(mean = c(0.5,0.5),
                sigma = matrix(c(1,0.8,0.8,1), nrow=2)) * 20   # multiply for more slope!
    })
}

# n = 1000:
# true sd: 1.66...
#  - without X variation: 1.41...
# wild bootstrap estimated: 0.606

# n = 500:
dgp <- function() {
  gen_df_kang(n = 500)
}


# true sd of atthat-att: 0.52
# wild bootstrap estimated: 0.24
dgp <- function(x=1) {
  gen_df_acic(
    model.trt="step", 
    root.trt=0.35, 
    overlap.trt="full",
    model.rsp="step", 
    alignment=0.75, 
    te.hetero="high",
    random.seed=x,        # NOTE: set random seed!
    n=1000, 
    p=10
  )
}

df_toy <- dgp()
dist_scaling_toy <- df_toy %>%
  summarize(across(starts_with("X"),
                   function(x) {
                     if (is.numeric(x)) 5 / (max(x) - min(x))
                     else 1000
                   }))

# record true sds ---------------------------------------------------------

# just do whole dgp again
#  - nice if true tau has zero variance, which it does here.

taus_resamp <- map_dbl(
  1:25,
  function(x) {
    d <- dgp() %>% 
      mutate(id = 1:n())
    
    att_hat <- d %>% 
      get_cal_matches(
        metric = "maximum",
        dist_scaling = dist_scaling_toy,
        cal_method = "fixed",
        est_method = "scm",
        return = "agg_co_units",
        knn = 25) %>% 
      get_att_ests()
    
    att <- d %>% 
      filter(Z) %>% 
      summarize(att = mean(Y1-Y0)) %>% 
      pull(att)
    
    print(paste("ATT:", round(att,3), "ATThat:", round(att_hat,3)))
    return(att_hat-att)
  }
)
sd(taus_resamp)


# conditional on simulated X, resample Z
taus_resampZ <- map_dbl(
  1:25,
  function(x) {
    d <- df_toy %>%    # put df_toy here! just resample its treatment
      rowwise() %>% 
      mutate(Z = sample(c(T,F), size=1, replace=T, 
                        prob=c(e,1-e))) %>% 
      ungroup() %>% 
      mutate(id = 1:n, .before=V1)
    
    att_hat <- d %>% 
      get_cal_matches(
        metric = "maximum",
        dist_scaling = dist_scaling_toy,
        cal_method = "fixed",
        est_method = "scm",
        return = "agg_co_units",
        knn = 25) %>% 
      get_att_ests()
    
    att <- 0
    print(paste("ATT:", round(att,3), "ATThat:", round(att_hat,3)))
    return(att_hat-att)
  }
)
sd(taus_resampZ)


# get bootstrapped estimates ----------------------------------------------

# get CSM results, aggregated co units
preds_csm <- get_cal_matches(
  df = df_toy,
  metric = "maximum",
  dist_scaling = dist_scaling_toy,
  cal_method = "fixed",
  est_method = "scm",
  return = "sc_units",
  knn = 25)   # NOTE: knn does nothing when cal_method = "fixed"...

boot_wild <- preds_csm %>% 
  wild_bootstrap(B=500)

sd(boot_wild)






# what if we just boot the whole thing?
boot_whole_thing_samps <- map_dbl(
  1:100,
  function(x) {
    d <- df_toy %>% 
      sample_n(size=n(), replace=T)
    
    d_scaling <- d %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) 6 / (max(x) - min(x))
                         else 1000
                       }))
    
    d %>% 
      get_cal_matches(
        metric = "maximum",
        dist_scaling = d_scaling,
        cal_method = "fixed",
        est_method = "scm",
        return = "agg_co_units",
        knn = 25) %>% 
      get_att_ests()
  }
)
sd(boot_whole_thing_samps)   # only 2x too large... just super slow...




# still super wrong.
if (F) {
  wf <- function(n) as.numeric(rmultinom(1, size=n, prob=rep(1/n,n)))
  # wf <- function(n) as.numeric(gtools::rdirichlet(1, alpha=rep(1,n)))
  
  wboot2_samps <- weighted_boot_naive(preds_csm, wf, B=500)
  wboot3_samps <- weighted_boot_norm(preds_csm, wf, B=500)

  sd(wboot2_samps)
  sd(wboot3_samps)
  sd(taus_resamp)
}



