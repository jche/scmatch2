
# attempting to write a simulation that breaks balancing,
#  a la toy example 2 in the paper

# idea: 
#  - control x-values mostly clustered on off-diagonal corners
#  - treatment x-values mostly clustered on diagonal corners
#  - a few control x-values near treatment x-values

require(tidyverse)
require(mvtnorm)



# simplest example in 2D, just like toy example
#  - constant treatment effect
gen_df_adv2d <- function(nc, nt, sig=0.2) {
  MU <- c(0.5, 0.5)
  SIG <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow=2)
  
  # tx units clustered at (0,0) and (1,1)
  dat_txblobs <- tibble(
    X1 = c(rnorm(nt, mean=0, sd=sig), rnorm(nt, mean=1, sd=sig)),
    X2 = c(rnorm(nt, mean=0, sd=sig), rnorm(nt, mean=1, sd=sig)),
    Z  = T
  )
  
  # co units clustered at (1,0) and (0,1)
  dat_coblobs <- tibble(
    X1 = c(rnorm(nc, mean=0, sd=sig), rnorm(nc, mean=1, sd=sig)),
    X2 = c(rnorm(nc, mean=1, sd=sig), rnorm(nc, mean=0, sd=sig)),
    Z  = F
  )
  
  # some co units near (0,0) and (1,1)
  dat_conear <- tibble(
    X1 = c(rnorm(nt, mean=0, sd=sig), rnorm(nt, mean=1, sd=sig)),
    X2 = c(rnorm(nt, mean=0, sd=sig), rnorm(nt, mean=1, sd=sig)),
    Z  = F
  )
  
  dat <- bind_rows(dat_txblobs, dat_coblobs, dat_conear)
  
  # function for control POs
  effect_fun <- function(X1, X2) { abs(X1-X2) }
  # constant tx effect
  tx_effect <- 0.2
  
  dat %>% 
    mutate(Y0 = effect_fun(X1,X2),
           Y1 = effect_fun(X1,X2) + tx_effect,
           Y  = ifelse(Z, Y1, Y0)) %>% 
    mutate(id = 1:n(), .before=X1)
  # rowwise() %>% 
  # mutate(Y0 = dmvnorm(x = c(X1,X2),
  #                     mean = MU,
  #                     sigma = SIG)) %>% 
  # ungroup()
}

if (F) {
  set.seed(90210)
  gen_df_adv2d(nc=100, nt=10, sig=0.2) %>% 
    ggplot(aes(X1,X2)) +
    geom_point(aes(pch=as.factor(Z), color=Y0))
}



# generate data -----------------------------------------------------------

# GOAL: predict Y0 for treated units

# generate data
set.seed(90210)
df <- gen_df_adv2d(nc=1000, nt=10, sig=0.2)
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


require(optweight)
m_bal <- optweight(Z ~ X1+X2,
                   data = df,
                   tols = c(0.01, 0.01),
                   estimand = "ATT")

# grab weights for each unit
preds_bal <- df %>% 
  mutate(wt = m_bal$weights)

# get ATT estimate: biased!
preds_bal %>% 
  group_by(Z) %>% 
  summarize(Y = mean(Y*wt))

# check balance: it's great!
if (F) {
  preds_bal %>% 
    group_by(Z) %>% 
    summarize(across(contains("X"),
                     ~mean(.x*wt)))
}


# BART --------------------------------------------------------------------

# History of packages (see Nonparametric Machine Learning and Efficient
#  Computation with Bayesian Additive Regression Trees: The BART R Package):
#  - BayesTree
#  - bartMachine
#  - dbarts ("drop-in replacement" for BayesTree, better with serial computation)
#  - BART (better with multi-threading)

# Hill uses BayesTree, so I'll use dbarts
# TODO: add a propensity score for even better results
require(dbarts)

m_bart <- bart(x.train = df %>% select(Z,X1,X2), 
               y.train = df$Y,
               x.test  = df %>% 
                 select(Z,X1,X2) %>% 
                 filter(Z==1) %>% 
                 mutate(Z=0) )

# grab counterfactual predictions for each treated unit
preds_bart <- df %>% 
  filter(Z==1) %>% 
  mutate(Y0hat = colMeans(m_bart$yhat.test))

# get ATT estimate: it's great!
preds_bart %>% 
  summarize(ATThat = mean(Y-Y0hat))

# visualize BART predictions
if (F) {
  # TODO
}



# BART propensity score ---------------------------------------------------

# TODO?


# our method --------------------------------------------------------------

source("R/distance.R")
source("R/sc.R")
source("R/matching.R")

preds_csm <- get_cal_matches(
  df = df,
  metric = "maximum",
  cal_method = "cem",
  est_method = "scm",
  return = "sc_units",
  num_bins = 5,
  wider = F,
  cem_tx_units = df %>% 
    filter(Z==1) %>% 
    pull(id)
)

# get ATT estimate:
preds_csm %>% 
  group_by(Z) %>% 
  summarize(Y = mean(Y*weights))

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
  method = "scm",
  return = "sc_units")

# get ATT estimate:
preds_cem %>% 
  group_by(Z) %>% 
  summarize(Y = mean(Y*weights))

# check balance: it's great!
if (F) {
  preds_cem %>% 
    group_by(Z) %>% 
    summarize(across(contains("X"),
                     ~mean(.x*weights)))
}





# TODO: some standard doubly robust thing? perhaps unnecessary
#  - see https://cran.r-project.org/web/packages/drtmle/vignettes/using_drtmle.html
#  - https://multithreaded.stitchfix.com/blog/2021/07/23/double-robust-estimator/



# adaptive hyperboxes -----------------------------------------------------

# TODO: some almost matching exactly thing?
#  - see https://almost-matching-exactly.github.io/software
#  - issue: MALTS is only in python, adaptive hyperboxes are in R though












