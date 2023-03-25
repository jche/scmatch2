
# Same as sim_adversarial, but it's actually the toy example used in the paper!

require(tidyverse)
require(mvtnorm)

require(optweight)
require(dbarts)
require(tmle)
require(AIPW)

require(tictoc)

source("R/distance.R")
source("R/sc.R")
source("R/matching.R")
source("R/estimate.R")
source("R/sim_data.R")
source("R/wrappers.R")
source("R/utils.R")


# set superlearner libraries
SL.library1 <- c("SL.mean", "SL.lm", "SL.glm")

# "based on Dorie et al. ACIC 2016 competition
#   and SuperLearner documentation,
#   adjusted to reduce the computational burden"
SL.library2 <- c("SL.glm", "SL.gam", "SL.glmnet",
                 # "SL.gbm", "SL.bartMachine",   # too slow
                 "SL.randomForest", "SL.xgboost")             # based on aipw recommendations

# tmle defaults
SL.library3Q <- c("SL.glm", "tmle.SL.dbarts2", "SL.glmnet")   # default tmle Q.SL.library
SL.library3g <-  c("SL.glm", "tmle.SL.dbarts.k.5", "SL.gam")  # default tmle g.SL.library



# run simulations ---------------------------------------------------------

zform1 <- as.formula("Z ~ X1+X2")
zform2 <- as.formula("Z ~ X1*X2")
form1 <- as.formula("Y ~ X1+X2")
form2 <- as.formula("Y ~ X1*X2")

# repeatedly run for all combinations of pars
tic()
for (i in 1:250) {
  nc <- 1000
  nt <- 200
  f0_sd <- 0.1
  df <- gen_df_adv(
    nc=nc, 
    nt=nt, 
    f0_sd = f0_sd,
    tx_effect_fun = function(X1, X2) {0.2},
    # f0_fun = function(x,y) {abs(x-y)})
    f0_fun = function(x,y) {
      matrix(c(x,y), ncol=2) %>%
        dmvnorm(mean = c(0.5,0.5),
                sigma = matrix(c(1,0.8,0.8,1), nrow=2)) * 4   # multiply for more slope!
      })
  
  dist_scaling <- df %>%
    summarize(across(starts_with("X"),
                     function(x) {
                       if (is.numeric(x)) 8 / (max(x) - min(x))
                       else 1000
                     }))

  res <- tibble(
    runid = i,
    nc = nc,
    nt = nt,
    f0_sd = f0_sd,

    true_ATT = df %>%
      filter(Z) %>%
      summarize(att = mean(Y1-Y0)) %>%
      pull(att),
    
    diff = get_att_diff(df),

    bal1 = get_att_bal(df, zform1, c(0.01, 0.01)),
    bal2 = get_att_bal(df, zform2, c(0.01, 0.01, 0.1)),

    or_lm = get_att_or_lm(df, form=form2),
    or_bart = get_att_or_bart(df, covs=c(X1,X2)),
    ps_lm = get_att_ps_lm(df, zform2),
    ps_bart = get_att_ps_bart(df, covs=c(X1,X2)),
    
    csm_scm = get_att_csm(df, dist_scaling=dist_scaling, est_method="scm"),
    csm_avg = get_att_csm(df, dist_scaling=dist_scaling, est_method="average"),
    cem_scm = get_att_cem(df, num_bins=8, est_method="scm"),
    cem_avg = get_att_cem(df, num_bins=8, est_method="average"),
    onenn = get_att_1nn(df, dist_scaling=dist_scaling),
    
    tmle1 = get_att_tmle(df %>% mutate(X3 = X1*X2), 
                         covs=c(X1,X2,X3),
                         Q.SL.library = SL.library1,
                         g.SL.library = SL.library1),
    # tmle2 = get_att_tmle(df, covs=c(X1,X2),
    #                      Q.SL.library = SL.library3Q,
    #                      g.SL.library = SL.library3g),
    # tmle3 = get_att_tmle(df, covs=c(X1,X2,X3),
    #                      Q.SL.library = SL.library3Q,
    #                      g.SL.library = SL.library3g),

    aipw1 = get_att_aipw(df %>% mutate(X3 = X1*X2), 
                         covs=c(X1,X2,X3),
                         Q.SL.library = SL.library1,
                         g.SL.library = SL.library1),
    # aipw2 = get_att_aipw(df, covs=c(X1,X2),
    #                      Q.SL.library = SL.library2,
    #                      g.SL.library = SL.library2),
    # aipw3 = get_att_aipw(df, covs=c(X1,X2),
    #                      Q.SL.library = SL.library2,
    #                      g.SL.library = SL.library2)
  )
  
  res %>% 
    pivot_longer(-c(runid:true_ATT))
  
  FNAME <- "sim_toy_results/toy_constant-tx-effect6.csv"
  if (file.exists(FNAME)) {
    write_csv(res, FNAME, append=T)
  } else {
    write_csv(res, FNAME)
  }
}
toc()


# analyze results ---------------------------------------------------------

# plot sample data
set.seed(1)
df <- gen_df_adv(
  nc=1000, 
  nt=100, 
  f0_sd = 0.1,
  tx_effect_fun = function(X1, X2) {0.2},
  # f0_fun = function(x,y) {abs(x-y)})
  f0_fun = function(x,y) {
    matrix(c(x,y), ncol=2) %>%
      dmvnorm(mean = c(0.5,0.5),
              sigma = matrix(c(1,0.8,0.8,1), nrow=2)) * 4   # multiply for more slope!
  })
vlines <- seq(min(df$X1),max(df$X1),length.out=9)
hlines <- seq(min(df$X2),max(df$X2),length.out=9)
df %>% 
  ggplot(aes(x=X1, y=X2, color=Y)) +
  geom_point(aes(pch=Z)) +
  scale_color_continuous(low="blue", high="orange") +
  geom_vline(xintercept=vlines, lty="dotted") +
  geom_hline(yintercept=hlines, lty="dotted") +
  facet_wrap(~Z)

if (F) {
  # for seed(1): cem does better! it drops a bunch though.
  df %>% 
    filter(Z) %>% 
    summarize(ATT = mean(Y1-Y0))
  csm_scm <- get_cal_matches(
    df = df,
    metric = "maximum",
    cal_method = "adaptive",
    est_method = "scm",
    dist_scaling = df %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) 8 / (max(x) - min(x))
                         else 1000
                       })),
    return = "all",
    knn = 25)
  csm_scm %>% 
    count(subclass) %>% 
    ggplot(aes(x=n)) +
    geom_histogram(binwidth=2, color="black")
  
  csm_avg = get_att_csm(df, num_bins=8, est_method="average")
  cem_scm = get_att_cem(df, num_bins=8, est_method="scm")
  cem_avg = get_att_cem(df, num_bins=8, est_method="average")
  
  csm_scm
  csm_avg
  cem_scm
  cem_avg
}




total_res <- read_csv("sim_toy_results/toy_constant-tx-effect6.csv")

total_res %>% 
  pivot_longer(-(runid:true_ATT)) %>%
  group_by(nc,nt,f0_sd,name) %>%
  summarize(rmse = sqrt(mean((value-true_ATT)^2)),
            bias = mean(value-true_ATT),
            sd   = sd(value)) %>% 
  rename(method = name) %>% 
  pivot_longer(c(rmse, bias, sd)) %>%
  ggplot(aes(x=method)) +
  geom_col(aes(y=value)) +
  facet_wrap(~name, scales="free") +
  coord_flip()

# prove identity
total_res %>% 
  pivot_longer(-(runid:true_ATT)) %>%
  group_by(nc,nt,f0_sd,name) %>%
  summarize(rmse = sqrt(mean((value-true_ATT)^2)),
            bias = mean(value-true_ATT),
            var_est = mean((value-mean(value))^2),
            var_true = mean((true_ATT-mean(true_ATT))^2),
            cov = mean((value-mean(value))*(true_ATT-mean(true_ATT)))) %>% 
  # mutate(rmse2 = sqrt(var_est + var_true + bias^2 - 2*cov),
  #        correct = all.equal(rmse,rmse2)) %>% 
  rename(method = name) %>% 
  pivot_longer(c(rmse, bias, var_est, var_true, cov)) %>% 
  ggplot(aes(x=method)) +
  geom_col(aes(y=value)) +
  facet_wrap(~name, scales="free") +
  coord_flip()

# results:
#  - toy_constant-tx-effect: 
#     - aipw1/bal2/csm_avg best...?
#     - 1nn lowest bias but very variable
#     - in general, does seem like lowering variance matters quite a bit...
#  - toy_constant-tx-effect2: doubled dnorm f0 function, for more slope there & cut noise
#     - 1nn and cem_scm suddenly very good, csm better than bal2/aipw1 but not great...
#     --> I think 1nn does very well with low noise.
#  - toy_constant-tx-effect3: triple dnorm, noise back up to 0.1
#     - still, 1nn and cem-scm doing very well in terms of bias for whatever reason...
#     - but csm better than everything else! moving in the right direction...?
#  - toy_constant-tx-effects4: quadruple dnorm, noise up to 0.2, smaller cells, more tx units
#     - pretty similar to 3...

# Q: HOW in the world does CEM outperform CSM in terms of bias?

#  - toy_constant-tx-effects5: same as 4, but fixed bug where 
#    csm is calipering on X3 as well (annoying that this changes things),
#    also added X3 to tmle estimator
#     - results: csm better now in terms of rmse, tmle1 not horrible anymore.
#     - challenges: avg is better than scm in terms of bias?? 
#                   cem is better than csm in terms of bias??
#  - toy_constant-tx-effects6: reduce noise again?
#     - doesn't seem to change results








