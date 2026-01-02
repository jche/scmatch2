
# attempting to write a simulation that breaks balancing,
#  a la toy example 2 in the paper

# idea: 
#  - control x-values mostly clustered on off-diagonal corners
#  - treatment x-values mostly clustered on diagonal corners
#  - a few control x-values near treatment x-values

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

tic()

#superlearner libraries

SL.library1 <- c("SL.mean", "SL.lm", "SL.glm")
# "based on Dorie et al. ACIC 2016 competition 
#   and SuperLearner documentation,
#   adjusted to reduce the computational burden"
SL.library2 <- c("SL.glm", "SL.gam", "SL.glmnet",
                 # "SL.gbm", "SL.bartMachine",   # too slow
                 "SL.randomForest", "SL.xgboost")




# simplest example in 2D, just like toy example
gen_df_adv2d <- function(nc, nt, 
                         tx_effect = 0.3,   # constant tx effect
                         sd = 0.1,
                         effect_fun = function(X1, X2) { abs(X1-X2) }) {
  
  # tx units clustered at (0,0) and (1,1)
  dat_txblobs <- tibble(
    X1 = c(rnorm(nt/2, mean=0, sd=0.2), 
           rnorm(nt/2, mean=1, sd=0.2)),
    X2 = c(rnorm(nt/2, mean=0, sd=0.2), 
           rnorm(nt/2, mean=1, sd=0.2)),
    Z  = T
  )
  
  # co units clustered at (1,0) and (0,1)
  dat_coblobs <- tibble(
    X1 = c(rnorm((nc-nt)/2, mean=0, sd=0.2), 
           rnorm((nc-nt)/2, mean=1, sd=0.2)),
    X2 = c(rnorm((nc-nt)/2, mean=1, sd=0.2), 
           rnorm((nc-nt)/2, mean=0, sd=0.2)),
    Z  = F
  )
  
  # some co units near (0,0) and (1,1)
  dat_conear <- tibble(
    X1 = c(rnorm(nt/2, mean=0, sd=0.2), 
           rnorm(nt/2, mean=1, sd=0.2)),
    X2 = c(rnorm(nt/2, mean=0, sd=0.2), 
           rnorm(nt/2, mean=1, sd=0.2)),
    Z  = F
  )
  
  dat <- bind_rows(dat_txblobs, dat_coblobs, dat_conear)
  
  res <- dat %>% 
    mutate(Y0 = effect_fun(X1,X2) + rnorm(n(), mean=0, sd=sd),
           Y1 = effect_fun(X1,X2) + tx_effect + rnorm(n(), mean=0, sd=sd),
           Y  = ifelse(Z, Y1, Y0)) %>% 
    mutate(id = 1:n(), .before=X1)
  print(paste0("SD of control outcomes: ", round(sd(res$Y0),4)))
  
  return(res)
}

if (F) {
  set.seed(90210)
  gen_df_adv2d(nc=100, nt=10, sig=0.2) %>% 
    ggplot(aes(X1,X2)) +
    geom_point(aes(pch=as.factor(Z), color=Y0))
}



# generate data -----------------------------------------------------------

# set.seed(90210)
df <- gen_df_adv2d(nc=1000, nt=10, 
                   tx_effect=0.2,
                   sd = 0.02,
                   effect_fun = function(x,y) {
                     matrix(c(x,y), ncol=2) %>% 
                       dmvnorm(mean = c(0.5,0.5),
                               sigma = matrix(c(1,0.8,0.8,1), nrow=2))
                   })
if (F) {
  df %>% 
    ggplot(aes(X1,X2)) +
    geom_point(aes(pch=as.factor(Z), color=Y0)) +
    scale_color_continuous(low="orange", high="blue") +
    theme_classic() +
    labs(pch = "Treated",
         color = latex2exp::TeX("$f_0(X)$"),
         x = latex2exp::TeX("$X_1$"),
         y = latex2exp::TeX("$X_2$"))
}



# balancing ---------------------------------------------------------------

# jose's package: sbw

# WeightIt package: https://rdrr.io/cran/WeightIt/man/weightit.html
#  - seems like method="optweight" gets us sbw?
#  - documentation suggests just using optweight package
# m_bal <- weightit(Z ~ X1+X2,
#                   data = df,
#                   estimand = "ATT",
#                   method = "optweight")


m_bal <- optweight(Z ~ X1+X2,
                   data = df,
                   tols = c(0.01, 0.01),
                   estimand = "ATT")

# grab weights for each unit
preds_bal <- df %>% 
  mutate(wt = m_bal$weights)

# get ATT estimate: biased!
ATT_bal <- preds_bal %>% 
  group_by(Z) %>% 
  summarize(Y = mean(Y*wt)) %>% 
  pivot_wider(names_from=Z, values_from=Y) %>% 
  mutate(ATThat = `TRUE` - `FALSE`) %>% 
  pull(ATThat)

# check balance: it's great!
if (F) {
  preds_bal %>% 
    group_by(Z) %>% 
    summarize(across(contains("X"),
                     ~mean(.x*wt)))
}


# outcome regression ------------------------------------------------------

# \frac{1}{n1} \sum Y_t - mhat0(X_t)

# See https://github.com/uber/causalml/issues/323 
#   for discussion on S vs. T learners
#  - S-learner: fit on full dataset
#  - T-learner: separately fit on just the co units
# --> I'll just use T-learner, doesn't drop much data (few tx units)

### lm

m_lm <- lm(Y ~ X1*X2, data = df %>% filter(!Z))
ATT_or_lm <- df %>% 
  filter(Z) %>% 
  mutate(mhat0 = predict(m_lm, newdata=.)) %>% 
  summarize(ATThat = mean(Y - mhat0)) %>% 
  pull(ATThat)


### BART

# History of packages (see Nonparametric Machine Learning and Efficient
#  Computation with Bayesian Additive Regression Trees: The BART R Package):
#  - BayesTree
#  - bartMachine
#  - dbarts ("drop-in replacement" for BayesTree, better with serial computation)
#  - BART (better with multi-threading)

# Hill uses BayesTree, so I'll use dbarts
# TODO: add a propensity score for even better results

m_bart <- bart(x.train = df %>% 
                 filter(Z==0) %>% 
                 select(X1,X2), 
               y.train = df %>% 
                 filter(Z==0) %>% 
                 pull(Y),
               x.test = df %>% 
                 filter(Z==1) %>% 
                 select(X1,X2))

# grab counterfactual predictions for each treated unit
ATT_or_bart <- df %>% 
  filter(Z==1) %>% 
  mutate(mhat0 = colMeans(m_bart$yhat.test)) %>% 
  summarize(ATThat = mean(Y-mhat0)) %>% 
  pull(ATThat)



# propensity score --------------------------------------------------------

# \frac{1}{n} \sum Y_t - e(Xj)/(1-e(Xj)) Y_j

### lm

m_lm_ps <- glm(Z ~ X1*X2, data = df, family="binomial")
ATT_ps_lm <- df %>% 
  mutate(e = predict(m_lm_ps, newdata=.),
         wt = ifelse(Z, 1, e/(1-e))) %>% 
  summarize(ATThat = sum(Z*wt*Y - (1-Z)*wt*Y) / n()) %>% 
  pull(ATThat)



# tmle --------------------------------------------------------------------

tmle1 <- tmle(Y = df$Y, 
              A = as.numeric(df$Z),
              W = df %>% 
                select(X1,X2) %>% 
                mutate(X3=X1*X2),
              Q.SL.library = SL.library1,
              g.SL.library = SL.library1)
ATT_tmle1 <- tmle1$estimates$ATT$psi

if (T) {
  tmle2 <- tmle(Y = df$Y, 
                A = as.numeric(df$Z),
                W = df %>% 
                  select(X1,X2) %>% 
                  mutate(X3=X1*X2),
                Q.SL.library = SL.library2,
                g.SL.library = SL.library2)
  ATT_tmle2 <- tmle2$estimates$ATT$psi
}


# timing
if (F) {
  require(tictoc)
  SL.library2 <- c("SL.glm", "SL.gam")
  tic()
  tmle2 <- tmle(Y = df$Y, 
                A = as.numeric(df$Z),
                W = df %>% 
                  select(X1,X2) %>% 
                  mutate(X3=X1*X2),
                Q.SL.library = SL.library2,
                g.SL.library = SL.library2)
  tmle2$estimates$ATT$psi
  toc()
  
  # glm: 0.59 sec
  # gam: 5.45 sec
  # glmnet: 7.5 sec
  # randomForest: 29.28 sec
  # xgboost: 64.51 sec
  # gbm: super slow, at least 5 mins
  # bartMachine: also super slow, at least 5 mins
}


# AIPW --------------------------------------------------------------------

# refresher on doubly robust estimation:
# http://www2.stat.duke.edu/~fl35/teaching/640/Chap3.5_Doubly%20Robust%20Estimation.pdf

# https://cran.r-project.org/web/packages/AIPW/AIPW.pdf
# https://github.com/yqzhong7/AIPW for short vignette

aipw1 <- AIPW$
  new(Y = df$Y,
      A = df$Z,
      W = df %>% 
        select(X1,X2) %>% 
        mutate(X3=X1*X2),
      Q.SL.library = SL.library1,
      g.SL.library = SL.library1,
      k_split = 5,
      verbose = T)$
  stratified_fit()$
  summary()
ATT_aipw1 <- aipw1$ATT_estimates$RD['Estimate'] %>% as.numeric()

# listWrappers()
# print(aipw1$result)
# aipw1$obs_est %>% names()

# takes a few minutes to run, is very accurate
if (T) {
  aipw2 <- AIPW$
    new(Y = df$Y,
        A = df$Z,
        W = df %>% select(X1,X2) %>% mutate(X3=X1*X2),
        Q.SL.library = SL.library2,
        g.SL.library = SL.library2,
        k_split = 5,
        verbose = T)$
    stratified_fit()$
    summary()
  ATT_aipw2 <- aipw2$ATT_estimates$RD['Estimate'] %>% as.numeric()
}


# timing
if (F) {
  require(tictoc)
  SL.library2 <- c("SL.glm", "SL.gam", "SL.glmnet", "SL.xgboost")
  tic()
  aipw2 <- AIPW$
    new(Y = df$Y,
        A = df$Z,
        W = df %>% select(X1,X2) %>% mutate(X3=X1*X2),
        Q.SL.library = SL.library2,
        g.SL.library = SL.library2,
        k_split = 5,
        verbose = T)$
    stratified_fit()$
    summary()
  aipw2$ATT_estimates$RD['Estimate']
  toc()
  
  # glm: 1.82 sec
  # gam: 12.5 sec, ish
  # glmnet: 19.2 sec, ish
  
  # randomForest: 134.35 sec
  # xgboost: 215 sec
}
  

# our method --------------------------------------------------------------

preds_csm <- get_cal_matches(
  df = df,
  metric = "maximum",
  cal_method = "fixed",
  est_method = "scm",
  dist_scaling = df %>%
    summarize(across(starts_with("X"),
                     function(x) {
                       if (is.numeric(x)) 5 / (max(x) - min(x))
                       else 1000
                     })),
  return = "sc_units",
  cem_tx_units = df %>% 
    filter(Z==1) %>% 
    pull(id)
)

# get ATT estimate:
ATT_csm <- preds_csm %>% 
  group_by(Z) %>% 
  summarize(Y = mean(Y*weights)) %>% 
  pivot_wider(names_from=Z, values_from=Y) %>% 
  mutate(ATThat = `TRUE` - `FALSE`) %>% 
  pull(ATThat)

# check balance: it's great!
if (F) {
  preds_csm %>% 
    group_by(Z) %>% 
    summarize(across(contains("X"),
                     ~mean(.x*weights)))
}


# CEM ---------------------------------------------------------------------

preds_cem <- get_cem_matches(
  df = df,
  num_bins = 5,
  est_method = "scm",
  return = "sc_units")

# get ATT estimate:
ATT_cem <- preds_cem %>% 
  group_by(Z) %>% 
  summarize(Y = mean(Y*weights)) %>% 
  pivot_wider(names_from=Z, values_from=Y) %>% 
  mutate(ATThat = `TRUE` - `FALSE`) %>% 
  pull(ATThat)

# check balance: it's great!
if (F) {
  preds_cem %>% 
    group_by(Z) %>% 
    summarize(across(contains("X"),
                     ~mean(.x*weights)))
}





# aggregate results -------------------------------------------------------

toc()   # ~12-13 minutes, ish

res <- c(
  true = 0.2,
  or_lm = ATT_or_lm,
  or_bart = ATT_or_bart,
  ps_lm = ATT_ps_lm,
  bal = ATT_bal,
  tmle1 = ATT_tmle1,
  tmle2 = ATT_tmle2,
  aipw1 = ATT_aipw1,
  aipw2 = ATT_aipw2,
  csm = ATT_csm,
  cem = ATT_cem)



# reproducing AIPW from TMLE ----------------------------------------------

# RESULTS: estimates are slightly different from TMLE
#  - part of this is random noise
#  - part of this is cross-fitting, I think...
# --> I'll just use AIPW instead of backing out from TMLE

if (F) {
  # AIPW ATT estimate 1
  df %>% 
    mutate(ehat = tmle1$g$g1W,
           mhat0 = as.numeric(tmle1$Qstar[,'Q0W'])) %>% 
    summarize(ATThat = 
                sum(Y*Z - 
                      (Y*(1-Z)*ehat + mhat0*(Z-ehat)) / (1-ehat) ) / 
                sum(Z) )
  
  # with initial Q estimate, not targeted one.
  df %>% 
    mutate(ehat = tmle1$g$g1W,
           mhat0 = as.numeric(tmle1$Qinit$Q[,'Q0W'])) %>% 
    summarize(ATThat = 
                sum(Y*Z - 
                      (Y*(1-Z)*ehat + mhat0*(Z-ehat)) / (1-ehat) ) / 
                sum(Z) )
  
  
  # AIPW ATT estimate 2
  df %>% 
    mutate(ehat = tmle1$g$g1W,
           mhat0 = as.numeric(tmle1$Qstar[,'Q0W'])) %>% 
    summarize(ATThat = 
                sum(Y*Z - 
                      (Y*(1-Z)*ehat + mhat0*(Z-ehat)) / (1-ehat) ) / 
                sum(Z) )
  
  
  
  
  df %>% 
    mutate(ehat = tmle1$g$g1W,
           mhat0 = as.numeric(tmle1$Qstar[,'Q0W'])) %>% 
    summarize(ATThat = 
                sum(Y*Z - 
                      (Y*(1-Z)*ehat + mhat0*(Z-ehat)) / (1-ehat) ) / 
                sum(Z) )
  
  
  # manually estimate ATE: 
  # using Qstar matches what AIPW package spits out
  df %>% 
    mutate(ehat = tmle1$g$g1W,
           mhat0 = as.numeric(tmle1$Qstar[,'Q0W']),
           mhat1 = as.numeric(tmle1$Qstar[,'Q1W'])) %>%    # with Qstar?
    summarize(ATEhat = 
                sum(Y*Z/ehat - (Z-ehat)*mhat1/ehat) / n() -
                sum(Y*(1-Z)/(1-ehat) + (Z-ehat)*mhat0/(1-ehat)) / n())
  
  # with AIPW package
  aipw_tmle1 <- AIPW_tmle$
    new(Y = df$Y, 
        A = as.numeric(df$Z),
        tmle_fit = tmle1,
        verbose = T)$
    summary(g.bound=0.025)
}



