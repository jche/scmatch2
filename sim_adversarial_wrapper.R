
# Idea: simplest possible simulation study where:
#  1) balancing performs poorly (if not balanced on enough things)
#  2) CSM performs on par with BART

# TODO: a more complex version of this study
#  - heterogeneous treatment effects?
#  - more covariates that interact like this?
#  - open Q: what would the goal of this more complex simulation be?
#     we've already shown (in this simulation) that there are cases 
#     where CSM is nice...

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


SL.library1 <- c("SL.mean", "SL.lm", "SL.glm")
# "based on Dorie et al. ACIC 2016 competition 
#   and SuperLearner documentation,
#   adjusted to reduce the computational burden"
SL.library2 <- c("SL.glm", "SL.gam", "SL.glmnet",
                 # "SL.gbm", "SL.bartMachine",   # too slow
                 "SL.randomForest", "SL.xgboost")


# functions for generating data -------------------------------------------

# simplest example in 2D, just like toy example
gen_df_adv2d <- function(nc, nt, 
                         tx_effect = 0.2,   # constant tx effect
                         effect_fun = function(X1, X2) { abs(X1-X2) },
                         sig = 0.2) {
  MU <- c(0.5, 0.5)
  SIG <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow=2)
  
  # tx units clustered at (0,0) and (1,1)
  dat_txblobs <- tibble(
    X1 = c(rnorm(nt/2, mean=0, sd=sig), 
           rnorm(nt/2, mean=1, sd=sig)),
    X2 = c(rnorm(nt/2, mean=0, sd=sig), 
           rnorm(nt/2, mean=1, sd=sig)),
    Z  = T
  )
  
  # co units clustered at (1,0) and (0,1)
  dat_coblobs <- tibble(
    X1 = c(rnorm((nc-nt)/2, mean=0, sd=sig), 
           rnorm((nc-nt)/2, mean=1, sd=sig)),
    X2 = c(rnorm((nc-nt)/2, mean=1, sd=sig), 
           rnorm((nc-nt)/2, mean=0, sd=sig)),
    Z  = F
  )
  
  # some co units near (0,0) and (1,1)
  dat_conear <- tibble(
    X1 = c(rnorm(nt/2, mean=0, sd=sig), 
           rnorm(nt/2, mean=1, sd=sig)),
    X2 = c(rnorm(nt/2, mean=0, sd=sig), 
           rnorm(nt/2, mean=1, sd=sig)),
    Z  = F
  )
  
  dat <- bind_rows(dat_txblobs, dat_coblobs, dat_conear)
  
  dat %>% 
    mutate(Y0 = effect_fun(X1,X2),
           Y1 = effect_fun(X1,X2) + tx_effect,
           Y  = ifelse(Z, Y1, Y0)) %>% 
    mutate(id = 1:n(), .before=X1)
}


# wrapper functions for getting att estimate ------------------------------

get_att_bal <- function(d, 
                        form,
                        tols) {
  m_bal <- optweight(form,
                     data = d,
                     tols = tols,
                     estimand = "ATT")
  
  # output ATT estimate
  d %>% 
    mutate(wt = m_bal$weights) %>% 
    group_by(Z) %>% 
    summarize(Y = mean(Y*wt)) %>% 
    pivot_wider(names_from=Z, names_prefix="Y", values_from=Y) %>% 
    mutate(ATT = YTRUE - YFALSE) %>% 
    pull(ATT)
}

get_att_or_lm <- function(d,
                          form) {
  m_lm <- lm(form, data = d %>% filter(!Z))
  
  d %>% 
    filter(Z) %>% 
    mutate(mhat0 = predict(m_lm, newdata=.)) %>% 
    summarize(ATThat = mean(Y - mhat0)) %>% 
    pull(ATThat)
}

get_att_or_bart <- function(d,
                         covs) {
  
  m_bart <- bart(x.train = d %>% 
                   filter(Z==0) %>% 
                   select({{covs}}), 
                 y.train = d %>% 
                   filter(Z==0) %>% 
                   pull(Y),
                 x.test = d %>% 
                   filter(Z==1) %>% 
                   select({{covs}}))
  
  # output ATT estimate
  d %>% 
    filter(Z==1) %>% 
    mutate(mhat0 = colMeans(m_bart$yhat.test)) %>% 
    summarize(ATThat = mean(Y-mhat0)) %>% 
    pull(ATThat)
}

get_att_ps_lm <- function(d,
                          form) {
  m_lm_ps <- glm(form, data = d, family="binomial")
  d %>% 
    mutate(e = predict(m_lm_ps, newdata=.),
           wt = ifelse(Z, 1, e/(1-e))) %>% 
    summarize(ATThat = sum(Z*wt*Y - (1-Z)*wt*Y) / n()) %>% 
    pull(ATThat)
}

get_att_tmle <- function(d, 
                         covs,
                         SL.library) {
  tmle <- tmle(Y = d$Y, 
               A = as.numeric(d$Z),
               W = d %>% 
                 select({{covs}}),
               Q.SL.library = SL.library,
               g.SL.library = SL.library)
  tmle$estimates$ATT$psi
}

get_att_aipw <- function(d,
                         covs,
                         SL.library) {
  aipw <- AIPW$
    new(Y = d$Y,
        A = d$Z,
        W = d %>% 
          select({{covs}}),
        Q.SL.library = SL.library1,
        g.SL.library = SL.library1,
        k_split = 5,
        verbose = F)$
    stratified_fit()$
    summary()
  
  aipw$ATT_estimates$RD['Estimate'] %>% 
    as.numeric()
}

get_att_csm <- function(d,
                        num_bins) {
  preds_csm <- get_cal_matches(
    df = d,
    metric = "maximum",
    cal_method = "cem",
    est_method = "scm",
    return = "sc_units",
    num_bins = num_bins,
    wider = F,
    cem_tx_units = d %>% 
      filter(Z==1) %>% 
      pull(id)
  )
  
  # get ATT estimate:
  preds_csm %>% 
    group_by(Z) %>% 
    summarize(Y = mean(Y*weights)) %>% 
    pivot_wider(names_from=Z, names_prefix="Y", values_from=Y) %>% 
    mutate(ATT = YTRUE - YFALSE) %>% 
    pull(ATT)
}

get_att_cem <- function(d,
                        num_bins) {
  preds_cem <- get_cem_matches(
    df = d,
    num_bins = num_bins,
    method = "average",
    return = "sc_units")
  
  # get ATT estimate:
  preds_cem %>% 
    group_by(Z) %>% 
    summarize(Y = mean(Y*weights)) %>% 
    pivot_wider(names_from=Z, names_prefix="Y", values_from=Y) %>% 
    mutate(ATT = YTRUE - YFALSE) %>% 
    pull(ATT)
}


# test these
if (F) {
  df <- gen_df_adv2d(nc=1000, nt=10, sig=0.2)
  df <- gen_df_adv2d(nc=1000, nt=10, 
                     tx_effect=0.2,
                     effect_fun = function(x,y) {
                       matrix(c(x,y), ncol=2) %>% 
                         dmvnorm(mean = c(0.5,0.5),
                                 sigma = matrix(c(1,0.8,0.8,1), nrow=2))
                     },
                     sig=0.2)
  
  get_att_bal(df, 
              form=as.formula('Z ~ X1+X2'), 
              tols=c(0.01, 0.01))
  get_att_bal(df, 
              form=as.formula('Z ~ X1*X2'), 
              tols=c(0.01, 0.01, 0.1))
  
  get_att_or_lm(df, form=as.formula('Y ~ X1*X2'))
  get_att_or_bart(df, covs=c(X1, X2))
  get_att_ps_lm(df, form=as.formula('Z ~ X1*X2'))
  
  get_att_tmle(df %>% 
                 mutate(X3 = X1*X2), 
               covs=c(X1,X2,X3),
               SL.library = SL.library1)
  get_att_aipw(df %>% 
                 mutate(X3 = X1*X2), 
               covs=c(X1,X2,X3),
               SL.library = SL.library1)
  
  get_att_csm(df, num_bins=5)
  get_att_cem(df, num_bins=5)
}



# run simulations ---------------------------------------------------------

pars <- expand_grid(
  nc = 1000,
  nt = 10,
  sig = 0.2
)

zform1 <- as.formula("Z ~ X1+X2")
zform2 <- as.formula("Z ~ X1*X2")
form1 <- as.formula("Y ~ X1+X2")
form2 <- as.formula("Y ~ X1*X2")

# repeatedly run for all combinations of pars
tic()
total_res <- map_dfr(
  1:3,
  function(x) {
    pars %>% 
      mutate(runid = x, .before=nc) %>% 
      rowwise() %>% 
      mutate(df = list(gen_df_adv2d(nc, nt, sig,
                                    effect_fun = function(x,y) {
                                      matrix(c(x,y), ncol=2) %>% 
                                        dmvnorm(mean = c(0.5,0.5),
                                                sigma = matrix(c(1,0.8,0.8,1), nrow=2))
                                    }) %>% 
                         mutate(X3 = X1*X2))) %>% 
      mutate(bal1 = get_att_bal(df, zform1, c(0.01, 0.01)),
             bal2 = get_att_bal(df, zform2, c(0.01, 0.01, 0.01)),
             
             or_lm = get_att_or_lm(df, form=form2),
             or_bart = get_att_or_bart(df, covs=c(X1, X2)),
             ps_lm = get_att_ps_lm(df, zform2),
             
             tmle1 = get_att_tmle(df, covs=c(X1,X2,X3), SL.library1),
             tmle2 = get_att_tmle(df, covs=c(X1,X2), SL.library2),
             tmle3 = get_att_tmle(df, covs=c(X1,X2,X3), SL.library2),
             
             aipw1 = get_att_aipw(df, covs=c(X1,X2,X3), SL.library1),
             aipw2 = get_att_aipw(df, covs=c(X1,X2), SL.library2),
             aipw3 = get_att_aipw(df, covs=c(X1,X2,X3), SL.library2),
             
             csm = get_att_csm(df, num_bins=5),
             cem = get_att_cem(df, num_bins=5)) %>% 
      ungroup() %>% 
      select(-df)
  })
toc()

write_csv(total_res, "sim_adversarial_results/test2.csv")




# analyze results ---------------------------------------------------------

# plot sample data
gen_df_adv2d(nc=1000, nt=10, sig=0.2) %>% 
  ggplot(aes(x=X1, y=X2, color=Y)) +
  geom_point(aes(pch=Z)) +
  scale_color_continuous(low="blue", high="orange")


total_res <- read_csv("sim_adversarial_results/test2.csv")
# show results
total_res %>% 
  pivot_longer(bal1:cem) %>% 
  group_by(name) %>% 
  summarize(rmse = sqrt(mean((value-0.2)^2)),
            bias = mean(value-0.2),
            sd   = sd(value))




