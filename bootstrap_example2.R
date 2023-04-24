
# checking what exactly starts to break the bootstrap
# - conclusion: increasing variance in f0(x) breaks the bootstrap
# --> this is because it's designed to average elements with mean zero,
#     not elements with mean f0(x)!


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

# wboots breaks for this (10x+ true sd of 0.1?)
#  - still happens if tx_effect_fun is 1
#  - much less bad if f0_fun isn't multiplied by 20 (only 2-3x true sd of 0.1)
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

# wboots not great for this (3x true sd of 0.1)
#  - idea: it's the f0_fun that's breaking things!
dgp <- function() {
  gen_df_adv(
    nc=500,
    nt=100,
    f0_sd = 0.5,
    tx_effect_fun = function(X1, X2) {1},
    f0_fun = function(x,y) {x + y})
}

# wboots pretty good for this (roughly equal to true sd of 0.1)
dgp <- function() {
  gen_df_adv(
    nc=500,
    nt=100,
    f0_sd = 0.5,
    tx_effect_fun = function(X1, X2) {1},
    f0_fun = function(x,y) {0})
}

# wboots pretty good (roughly equal to true sd of 0.1)
dgp <- function() {
  gen_df_adv(
    nc=500,
    nt=100,
    f0_sd = 0.5,
    tx_effect_fun = function(X1, X2) {X1+X2},
    f0_fun = function(x,y) {0})
}


# # wboots work for this: constant tx effect, some imbalance
# dgp <- function() {
#   numsamp <- 500
#   df_toy <- tibble(
#     X1 = rnorm(numsamp),
#     X2 = rnorm(numsamp),
#     tau = 1,
#     e = invlogit(X1+X2-2),
#     Z = rbernoulli(n=numsamp, p=e),
#     Y = tau*Z + rnorm(numsamp)
#   ) %>% 
#     mutate(id = 1:n())
# }

# # wboots work for this too: roughly true sd
# dgp <- function() {
#   numsamp <- 500
#   df_toy <- tibble(
#     X1 = rnorm(numsamp),
#     X2 = rnorm(numsamp),
#     tau = X1+X2,
#     e = invlogit(X1+X2-2),
#     Z = rbernoulli(n=numsamp, p=e),
#     Y = tau*Z + rnorm(numsamp)
#   ) %>%
#     mutate(id = 1:n())
# }

df_toy <- dgp()
dist_scaling_toy <- df_toy %>%
  summarize(across(starts_with("X"),
                   function(x) {
                     if (is.numeric(x)) 6 / (max(x) - min(x))
                     else 1000
                   }))

ggplot(df_toy, aes(X1,X2)) +
  geom_point(aes(color=Z))

# record true sds ---------------------------------------------------------

# just do whole dgp again
#  - nice if true tau has zero variance, which it does here.

taus_resamp <- map_dbl(
  1:25,
  function(x) {
    dgp() %>% 
      mutate(id = 1:n()) %>% 
      get_cal_matches(
        metric = "maximum",
        dist_scaling = dist_scaling_toy,
        cal_method = "fixed",
        est_method = "scm",
        return = "agg_co_units",
        knn = 25) %>% 
      get_att_ests()
  }
)
sd(taus_resamp)

# get bootstrapped estimates ----------------------------------------------

# get CSM results, aggregated co units
preds_csm <- get_cal_matches(
  df = df_toy,
  metric = "maximum",
  dist_scaling = dist_scaling_toy,
  cal_method = "fixed",
  est_method = "scm",
  return = "agg_co_units",
  knn = 25)   # NOTE: knn does nothing when cal_method = "fixed"...

if (F) {
  ggplot(preds_csm, aes(X1,X2)) +
    geom_point(aes(color=Z, size=weights))
  co_weights <- preds_csm %>% 
    filter(!Z) %>% 
    pull(weights)
  mean(co_weights)
  sd(co_weights)
  hist(co_weights)
}

wf <- function(n) as.numeric(rmultinom(1, size=n, prob=rep(1/n,n)))
# wf <- function(n) as.numeric(gtools::rdirichlet(1, alpha=rep(1,n)))
# wf <- function(n) {
#   sample(
#     c( -(sqrt(5)-1)/2, (sqrt(5)+1)/2 ),
#     prob = c( (sqrt(5)+1)/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5)) ),
#     replace = T, size = n)
# }
boot_fsatt_naive <- preds_csm %>%
  weighted_boot_naive(wt_fun = wf, B=500)
boot_fsatt_norm <- preds_csm %>%
  weighted_boot_norm(wt_fun = wf, B=500)
boot_wild <- preds_csm %>% 
  wild_bootstrap(B=500)

# hist(boot_fsatt_naive)
# hist(boot_fsatt_norm)
sd(boot_fsatt_naive)
sd(boot_fsatt_norm)
sd(boot_wild)






# what if we just boot the whole thing?
boot_whole_thing_samps <- map_dbl(
  1:25,
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










# ground truth ------------------------------------------------------------

# true sd: randomness only in Z
#  - weird w.r.t. the ATT in some sense...
B <- 500
taus_resampZ <- map_dbl(
  1:B,
  function(x) {
    df %>% 
      mutate(Z = sample(0:1, size=numsamp, replace=T, prob=c(0.9,0.1)),
             Yobs = ifelse(Z, Y1, Y0),
             wt = ifelse(Z, 1, sum(Z)/(n()-sum(Z)))) %>% 
      summarize(att = sum((2*Z-1) * wt * Yobs) / sum(Z)) %>% 
      pull(att)
  })
ggplot(tibble(x=taus_resampZ)) +
  geom_histogram(aes(x)) +
  geom_vline(xintercept=tau)
mean(taus_resampZ)
sd(taus_resampZ)

# true sd: randomness in Y
taus_resampY <- map_dbl(
  1:B,
  function(x) {
    df %>% 
      mutate(Y0 = rnorm(numsamp),
             tauj = rnorm(numsamp, mean=tau),
             Y1 = Y0 + tauj,
             Yobs = ifelse(Z, Y1, Y0),
             wt = ifelse(Z, 1, sum(Z)/(n()-sum(Z)))) %>% 
      summarize(att = sum((2*Z-1) * wt * Yobs) / sum(Z)) %>% 
      pull(att)
  })
ggplot(tibble(x=taus_resampY)) +
  geom_histogram(aes(x)) +
  geom_vline(xintercept=tau)
mean(taus_resampY)
sd(taus_resampY)

# note: given original dataset may not be similar to the "average" dataset!



# TODO: what is the "true" variation we're after?
#  - in some sense, we'd like to condition on the treated units...?
# regardless, bootstrap seems to be much closer to resampZ se than to resampY se...




# bootstrap results -------------------------------------------------------

# TODO CHECK: does it matter what the co units' weights sum to...?

boot <- function(d, B=100) {
  map_dbl(1:B,
          .progress = "Bootstrapping...",
          function(b) {
            d %>% 
              sample_n(n(), replace=T) %>%
              summarize(att = sum((2*Z-1) * wt * Yobs) / sum(Z)) %>%
              pull(att)
          })
}

boot2 <- function(d, B=100) {
  # n <- nrow(d)
  # n1 <- sum(d$Z)
  
  map_dbl(1:B,
          .progress = "Bootstrapping...",
          function(b) {
            att_hat <- d %>%
              summarize(att = sum((2*Z-1) * wt * Yobs) / sum(Z)) %>%
              pull(att)
            
            d %>% 
              sample_n(n(), replace=T) %>%
              summarize(att = sum((2*Z-1) * wt * Yobs - sum(Z)/n()*att_hat) / sum(Z)) %>%
              # summarize(att = sum((2*Z-1) * wt * Yobs - n1/n*att_hat) / n1) %>%
              pull(att)
          })
}

# use n1/n*att_hat to cancel bias for multinomial weights,
#  which sum to n instead of n1
#  - this actually helps the wild bootstrap, for whatever reason???
# NOTE: this is conservative, due to n1 weirdness
weighted_boot <- function(d, wt_fun, B=100) {
  n <- nrow(d)
  n1 <- sum(d$Z)
  
  att_hat <- d %>%
    summarize(att = sum((2*Z-1) * wt * Yobs) / n1) %>%
    pull(att)
  
  map_dbl(1:B,
          .progress = "Bootstrapping...",
          function(b) {
            d %>% 
              mutate(
                weights_boot = wt_fun(n),
                tau_i = (2*Z-1) * wt * Yobs) %>% 
              summarize(att = sum(weights_boot * 
                                    (tau_i - n1/n*att_hat)) / n1) %>%
              pull(att)
          })
}

# this is pretty close to what the bootstrap does!
weighted_boot2 <- function(d, wt_fun, B=100) {
  att_hat <- d %>%
    summarize(att = sum((2*Z-1) * wt * Yobs) / sum(Z)) %>%
    pull(att)
  
  map_dbl(1:B,
          .progress = "Bootstrapping...",
          function(b) {
            d %>% 
              mutate(
                weights_boot = wt_fun(n()),
                tau_i = (2*Z-1) * wt * Yobs) %>% 
              summarize(
                # n1 = sum(weights_boot*Z),
                # n = sum(weights_boot),
                att = sum(weights_boot * 
                            (tau_i - sum(weights_boot*Z)/n()*att_hat)) / 
                  sum(weights_boot*Z) ) %>%
              pull(att)
          })
}

# this is pretty close to what the bootstrap does!
weighted_boot3 <- function(d, wt_fun, B=100) {
  map_dbl(1:B,
          .progress = "Bootstrapping...",
          function(b) {
            d %>% 
              mutate(
                weights_boot = wt_fun(n()),
                tau_i = (2*Z-1) * wt * Yobs) %>% 
              summarize(
                att = sum(weights_boot * (tau_i)) / sum(weights_boot*Z) ) %>%
              pull(att)
          })
}

if (F) {
  # full bootstrap works great
  boot_samps <- boot(df, B=500)
  hist(boot_samps)
  sd(boot_samps)
  
  boot_samps2 <- boot2(df, B=500)
  hist(boot_samps2)
  sd(boot_samps2)
  
  # weighted bootstrap works...?
  
  wf <- function(n) as.numeric(rmultinom(1, size=n, prob=rep(1/n,n)))
  # wf <- function(n) as.numeric(gtools::rdirichlet(1, alpha=rep(1,n)))
  # wf <- function(n) {
  #   sample(
  #     c( -(sqrt(5)-1)/2, (sqrt(5)+1)/2 ),
  #     prob = c( (sqrt(5)+1)/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5)) ),
  #     replace = T, size = n)
  # }
  
  wboot_samps <- weighted_boot(df, wf, B=500)
  hist(wboot_samps)
  sd(wboot_samps)
  
  wboot2_samps <- weighted_boot2(df, wf, B=500)
  hist(wboot2_samps)
  sd(wboot2_samps)
  
  wboot3_samps <- weighted_boot3(df, wf, B=500)
  hist(wboot3_samps)
  sd(wboot3_samps)
  
  sd(boot_samps)
  sd(boot_samps2)
  sd(wboot_samps)
  sd(wboot2_samps)
  sd(wboot3_samps)
  sd(taus_resampY)
  sd(taus_resampZ)
}



