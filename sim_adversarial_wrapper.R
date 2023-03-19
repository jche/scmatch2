
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


logit <- function(x) {
  log(x/(1-x))
}
invlogit <- function(x) {
  exp(x) / (1+exp(x))
}



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


# functions for generating data -------------------------------------------

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

# TODO: do this directly with mvnorms,
#  otherwise you get a carved-out "cross" not a carved-out "circle"
gen_df_veryadv <- function(nc, nt, 
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
  dcustom <- function(x, mean, sd) {
    dnorm(x, mean=mean, sd=sd) - dnorm(x, mean=mean, sd=sd/2)/2
  }
  rcustom <- function(n, mean, sd) {
    res <- c()
    while (length(res) < n) {
      proposals <- rnorm(n, mean, 2*sd)
      accept <- runif(n) < dcustom(proposals,mean,sd) / dnorm(proposals,mean,2*sd)
      accepted <- proposals[which(accept)]
      res <- c(res, accepted)
    }
    return(res[1:n])
  }
  if (F) {
    tibble(x = seq(-5,5,0.01)) %>% 
      mutate(dens = dcustom(x, mean=0, sd=1),
             dens2 = dnorm(x, mean=0, sd=2)) %>% 
      pivot_longer(dens:dens2) %>% 
      ggplot(aes(x=x)) +
      geom_histogram(data=tibble(x=rcustom(1000, 0, 1)),
                     aes(y=..density..)) +
      geom_point(aes(y=value*2, color=name))
  }
  
  dat_conear <- tibble(
    X1 = c(rcustom(nt/2, mean=0, sd=sig), 
           rcustom(nt/2, mean=1, sd=sig)),
    X2 = c(rcustom(nt/2, mean=0, sd=sig), 
           rcustom(nt/2, mean=1, sd=sig)),
    Z  = F
  )
  
  dat <- bind_rows(dat_txblobs, dat_coblobs, dat_conear)
  
  dat %>% 
    mutate(Y0 = effect_fun(X1,X2),
           Y1 = effect_fun(X1,X2) + tx_effect,
           Y  = ifelse(Z, Y1, Y0)) %>% 
    mutate(id = 1:n(), .before=X1)
}



# simplest example in 2D, just like toy example
gen_df_advhet <- function(nc, nt, 
                          tx_effect = function(X1, X2) {(X1-0.5)^2+(X2-0.5)^2},
                          eps_sd = 0.1,
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
    mutate(Y0 = effect_fun(X1,X2) + rnorm(n(), mean=0, sd=eps_sd),
           Y1 = effect_fun(X1,X2) + tx_effect(X1,X2) + rnorm(n(), mean=0, sd=eps_sd),
           Y  = ifelse(Z, Y1, Y0)) %>% 
    mutate(id = 1:n(), .before=X1)
  print(paste0("SD of control outcomes: ", round(sd(res$Y0),4)))
  
  return(res)
}




# generate co in unit square
gen_df_full <- function(nc, nt, eps_sd = 0.1,
                        tx_effect = function(X1, X2) {(X1-0.5)^2+(X2-0.5)^2},
                        effect_fun = function(X1, X2) { abs(X1-X2) }) {
  dat_co <- tibble(
    X1 = runif(nc),
    X2 = runif(nc),
    Z = F
  )
  
  dat_tx1 <- tibble(
    X1 = rnorm(nt/2, mean=0.25, sd=0.1),
    X2 = rnorm(nt/2, mean=0.25, sd=0.1),
    Z = T
  )
  dat_tx2 <- tibble(
    X1 = rnorm(nt/2, mean=0.75, sd=0.1),
    X2 = rnorm(nt/2, mean=0.75, sd=0.1),
    Z = T
  )
  
  dat <- bind_rows(dat_co, dat_tx1, dat_tx2)
  print(dat %>% mutate(eff = tx_effect(X1,X2)) %>% head())
  res <- dat %>% 
    mutate(Y0 = effect_fun(X1,X2) + rnorm(n(), mean=0, sd=eps_sd),
           Y1 = effect_fun(X1,X2) + tx_effect(X1,X2) + rnorm(n(), mean=0, sd=eps_sd),
           Y  = ifelse(Z, Y1, Y0)) %>% 
    mutate(id = 1:n(), .before=X1)
  print(paste0("SD of control outcomes: ", round(sd(res$Y0),4)))
  
  return(res)
}


gen_df_ps <- function(n, eps_sd = 0.1,
                      ps = function(x,y) {invlogit(2*x+2*y - 5)},
                      tx_effect = function(X1, X2) {(X1-0.5)^2+(X2-0.5)^2},
                      effect_fun = function(X1, X2) { abs(X1-X2) }) {
  dat <- tibble(
    X1 = runif(n),
    X2 = runif(n)
  ) %>% 
    mutate(Z = runif(n) < ps(X1,X2))
  
  res <- dat %>% 
    mutate(Y0 = effect_fun(X1,X2) + rnorm(n(), mean=0, sd=eps_sd),
           Y1 = effect_fun(X1,X2) + tx_effect(X1,X2) + rnorm(n(), mean=0, sd=eps_sd),
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
  
  # TODO: WHY IS PS NOT WORKING?
  
  d %>% 
    mutate(e = invlogit(predict(m_lm_ps, newdata=.)),
           wt = ifelse(Z, 1, e/(1-e))) %>% 
    summarize(ATThat = sum(Z*wt*Y - (1-Z)*wt*Y) / n()) %>% 
    pull(ATThat)
}

get_att_ps_bart <- function(d,
                            covs) {
  
  m_bart <- bart(x.train = d %>% 
                   select({{covs}}), 
                 y.train = d$Z %>% as.numeric(),
                 x.test = d %>% 
                   select({{covs}}))
  
  # Q: why is bart not even getting ps right? shape is right, but values are low
  if (F) {
    gen_df_ps
    d %>% 
      mutate(e = pnorm(colMeans(m_bart$yhat.test))) %>% 
      ggplot(aes(X1,X2)) +
      geom_point(aes(color=e, shape=Z)) +
      facet_wrap(~Z)
  }
  
  # output ATT estimate
  d %>% 
    mutate(e = pnorm(colMeans(m_bart$yhat.test)),
           wt = Z - e*(1-Z)/(1-e)) %>% 
    summarize(ATThat = mean(wt*Y)) %>% 
    pull(ATThat)
}

get_att_tmle <- function(d, 
                         covs,
                         Q.SL.library,
                         g.SL.library) {
  tmle <- tmle(Y = d$Y, 
               A = as.numeric(d$Z),
               W = d %>% 
                 select({{covs}}),
               Q.SL.library = Q.SL.library,
               g.SL.library = g.SL.library,
               V = 5)
  tmle$estimates$ATT$psi
}

get_att_aipw <- function(d,
                         covs,
                         Q.SL.library,
                         g.SL.library) {
  aipw <- AIPW$
    new(Y = d$Y,
        A = d$Z,
        W = d %>% 
          select({{covs}}),
        Q.SL.library = Q.SL.library,
        g.SL.library = g.SL.library,
        k_split = 5,
        verbose = F)$
    stratified_fit()$
    summary()
  
  aipw$ATT_estimates$RD['Estimate'] %>% 
    as.numeric()
}

get_att_csm <- function(d,
                        num_bins,
                        est_method = "scm") {
  preds_csm <- get_cal_matches(
    df = d,
    metric = "maximum",
    cal_method = "adaptive",
    est_method = est_method,
    dist_scaling = d %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) num_bins / (max(x) - min(x))
                         else 1000
                       })),
    return = "sc_units",
    knn = 25)
  
  # get ATT estimate:
  preds_csm %>% 
    group_by(Z) %>% 
    summarize(Y = mean(Y*weights)) %>% 
    pivot_wider(names_from=Z, names_prefix="Y", values_from=Y) %>% 
    mutate(ATT = YTRUE - YFALSE) %>% 
    pull(ATT)
}

# estimates ATT!
get_att_cem <- function(d,
                        num_bins,
                        est_method = "average") {
  preds_feasible <- get_cem_matches(
    df = d,
    num_bins = num_bins,
    est_method = est_method,
    return = "sc_units")
  
  # get ATT estimate:
  att_feasible <- preds_feasible %>% 
    group_by(Z) %>% 
    summarize(Y = mean(Y*weights)) %>% 
    pivot_wider(names_from=Z, names_prefix="Y", values_from=Y) %>% 
    mutate(ATT = YTRUE - YFALSE) %>% 
    pull(ATT)
  
  if (length(attr(preds_feasible, "feasible_units")) < sum(d$Z)) {
    preds_infeasible <- d %>% 
      filter(!Z | !(id %in% attr(preds_feasible, "feasible_units"))) %>% 
      get_cal_matches(.,
                      metric = "maximum",
                      cal_method = "1nn",
                      est_method = est_method,
                      dist_scaling = d %>%
                        summarize(across(starts_with("X"),
                                         function(x) {
                                           if (is.numeric(x)) num_bins / (max(x) - min(x))
                                           else 1000
                                         })),
                      return = "sc_units")
    att_infeasible <- preds_infeasible %>% 
      group_by(Z) %>% 
      summarize(Y = mean(Y*weights)) %>% 
      pivot_wider(names_from=Z, names_prefix="Y", values_from=Y) %>% 
      mutate(ATT = YTRUE - YFALSE) %>% 
      pull(ATT)
    
    return((att_feasible * sum(preds_feasible$Z) + 
              att_infeasible * sum(preds_infeasible$Z)) / sum(d$Z))
  }
  
  return(att_feasible)
}


get_att_1nn <- function(d, num_bins) {
  preds_1nn <- get_cal_matches(
    df = d,
    metric = "maximum",
    cal_method = "1nn",
    est_method = "average",
    dist_scaling = d %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) num_bins / (max(x) - min(x))
                         else 1000
                       })),
    return = "sc_units",
    cem_tx_units = d %>% 
      filter(Z==1) %>% 
      pull(id)
  )
  
  # get ATT estimate:
  preds_1nn %>% 
    group_by(Z) %>% 
    summarize(Y = mean(Y*weights)) %>% 
    pivot_wider(names_from=Z, names_prefix="Y", values_from=Y) %>% 
    mutate(ATT = YTRUE - YFALSE) %>% 
    pull(ATT)
}



# test these --------------------------------------------------------------

if (F) {
  df <- gen_df_adv2d(nc=1000, nt=50, 
                     tx_effect=0.2,
                     sd=0.03,
                     effect_fun = function(x,y) {
                       matrix(c(x,y), ncol=2) %>% 
                         dmvnorm(mean = c(0.5,0.5),
                                 sigma = matrix(c(1,0.8,0.8,1), nrow=2))
                     })
  
  df <- gen_df_full(nc=1000, nt=50, eps_sd=0,
                    tx_effect = function(x,y) {(0.5*(x+y))},
                     effect_fun = function(x,y) {x+y})
  
  df <- gen_df_ps(n=10000, eps_sd=0, tx_effect = function(x,y) {1}, effect_fun = function(x,y) {1})
  
  df %>% 
    ggplot(aes(X1,X2)) +
    geom_point(aes(pch=as.factor(Z), color=Y)) +
    scale_color_continuous(low="orange", high="blue") +
    theme_classic() +
    labs(pch = "Treated",
         color = latex2exp::TeX("$f_0(X)$"),
         x = latex2exp::TeX("$X_1$"),
         y = latex2exp::TeX("$X_2$")) +
    facet_wrap(~Z)
  df %>% 
    filter(Z) %>% 
    summarize(eff = mean(Y1-Y0))
  
  
  get_att_bal(df,
              form=as.formula('Z ~ X1+X2'),
              tols=c(0.01, 0.01))
  get_att_bal(df,
              form=as.formula('Z ~ X1*X2'),
              tols=c(0.01, 0.01, 0.1))

  get_att_or_lm(df, form=as.formula('Y ~ X1*X2'))
  get_att_or_bart(df, covs=c(X1, X2))
  get_att_ps_lm(df, form=as.formula('Z ~ X1*X2'))
  get_att_ps_bart(df, covs=c(X1, X2))
  
  get_att_tmle(df %>% 
                 mutate(X3 = X1*X2), 
               covs=c(X1,X2,X3),
               Q.SL.library = SL.library1,
               g.SL.library = SL.library1)
  get_att_aipw(df %>% 
                 mutate(X3 = X1*X2), 
               covs=c(X1,X2,X3),
               Q.SL.library = SL.library1,
               g.SL.library = SL.library1)
  
  get_att_csm(df, num_bins=5, est_method="scm")
  get_att_cem(df, num_bins=5, est_method="scm")
}



# run simulations ---------------------------------------------------------

zform1 <- as.formula("Z ~ X1+X2")
zform2 <- as.formula("Z ~ X1*X2")
form1 <- as.formula("Y ~ X1+X2")
form2 <- as.formula("Y ~ X1*X2")

# repeatedly run for all combinations of pars
tic()
for (i in 1:250) {
  nc <- 1000
  nt <- 100
  eps_sd <- 0.05
  df <- gen_df_full(
    nc=nc, nt=nt, eps_sd = eps_sd,
    tx_effect = function(X1, X2) {(X1+X2)*2},
    effect_fun = function(x,y) {
      matrix(c(x,y), ncol=2) %>% 
        dmvnorm(mean = c(0.5,0.5),
                sigma = matrix(c(1,0.8,0.8,1), nrow=2))}) %>% 
    mutate(X3 = X1*X2)
  
  res <- tibble(
    runid = i,
    nc = nc,
    nt = nt,
    eps_sd = eps_sd,

    true_ATT = df %>%
      filter(Z) %>%
      summarize(att = mean(Y1-Y0)) %>%
      pull(att),

    bal1 = get_att_bal(df, zform1, c(0.01, 0.01)),
    bal2 = get_att_bal(df, zform2, c(0.01, 0.01, 0.1)),

    or_lm = get_att_or_lm(df, form=form2),
    or_bart = get_att_or_bart(df, covs=c(X1,X2)),
    ps_lm = get_att_ps_lm(df, zform2),
    ps_bart = get_att_ps_bart(df, covs=c(X1,X2)),
    
    csm_scm = get_att_csm(df, num_bins=5, est_method="scm"),
    csm_avg = get_att_csm(df, num_bins=5, est_method="average"),
    cem_scm = get_att_cem(df, num_bins=5, est_method="scm"),
    cem_avg = get_att_cem(df, num_bins=5, est_method="average"),
    onenn = get_att_1nn(df, num_bins=5),
    
    tmle1 = get_att_tmle(df, covs=c(X1,X2),
                         Q.SL.library = SL.library1,
                         g.SL.library = SL.library1),
    # tmle2 = get_att_tmle(df, covs=c(X1,X2),
    #                      Q.SL.library = SL.library3Q,
    #                      g.SL.library = SL.library3g),
    tmle3 = get_att_tmle(df, covs=c(X1,X2,X3),
                         Q.SL.library = SL.library3Q,
                         g.SL.library = SL.library3g),

    aipw1 = get_att_aipw(df, covs=c(X1,X2,X3),
                         Q.SL.library = SL.library1,
                         g.SL.library = SL.library1),
    # aipw2 = get_att_aipw(df, covs=c(X1,X2),
    #                      Q.SL.library = SL.library2,
    #                      g.SL.library = SL.library2),
    aipw3 = get_att_aipw(df, covs=c(X1,X2),
                         Q.SL.library = SL.library2,
                         g.SL.library = SL.library2)
  )
  
  FNAME <- "sim_adversarial_results/test_full4.csv"
  if (file.exists(FNAME)) {
    write_csv(res, FNAME, append=T)
  } else {
    write_csv(res, FNAME)
  }
}
toc()



# NOTES:
#  - test_full2: tx_effect = function(X1, X2) {(X1-0.5)^2+(X2-0.5)^2}
#  - test_full3: tx_effect = function(X1, X2) {X1+X2}
#  - test_full4: tx_effect = function(X1, X2) {(X1+X2)*2}, less sd
#  - 


# analyze results ---------------------------------------------------------

# plot sample data
set.seed(1)
# dat <- gen_df_advhet(nc=1000, nt=50)
dat <- gen_df_full(nc=1000, nt=100)
vlines <- seq(min(dat$X1),max(dat$X1),length.out=6)
hlines <- seq(min(dat$X2),max(dat$X2),length.out=6)
dat %>% 
  ggplot(aes(x=X1, y=X2, color=Y)) +
  geom_point(aes(pch=Z)) +
  scale_color_continuous(low="blue", high="orange") +
  geom_vline(xintercept=vlines, lty="dotted") +
  geom_hline(yintercept=hlines, lty="dotted")




total_res <- read_csv("sim_adversarial_results/test_full4.csv")
# show results

total_res %>% 
  pivot_longer(-runid) %>% 
  group_by(name) %>% 
  summarize(rmse = sqrt(mean((value-tx_effect)^2)),
            bias = mean(value-tx_effect),
            sd   = sd(value)) %>% 
  rename(method = name) %>% 
  # pivot_longer(c(rmse, bias, sd)) %>% 
  ggplot(aes(x=method)) +
  geom_col(aes(y=rmse)) +
  coord_flip()


total_res %>% 
  pivot_longer(-(runid:true_ATT)) %>%
  filter(!name %in% c("ps_lm", "ps_bart")) %>% 
  group_by(nc,nt,eps_sd,name) %>%
  summarize(rmse = sqrt(mean((value-true_ATT)^2)),
            bias = mean(value-true_ATT),
            sd   = sd(value)) %>% 
  rename(method = name) %>% 
  pivot_longer(c(rmse, bias, sd)) %>%
  ggplot(aes(x=method)) +
  geom_col(aes(y=value)) +
  facet_wrap(~name) +
  coord_flip()

# CEM does worse now, since it suffers some bias (due to dropping units with extreme tx effects)
# but 1nn still completely nails it! As does BART, but that's okay.

