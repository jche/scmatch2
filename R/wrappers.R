
# wrapper functions: input data/settings, output ATT estimate

require(optweight)
require(dbarts)
require(tmle)
require(AIPW)


get_att_diff <- function(d) {
  d %>% 
    group_by(Z) %>% 
    summarize(mn = mean(Y)) %>% 
    pivot_wider(names_from=Z, names_prefix="Y", values_from=mn) %>% 
    mutate(ATT = YTRUE - YFALSE) %>% 
    pull(ATT)
}

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
    mutate(e = invlogit(predict(m_lm_ps, newdata=.)),
           wt = ifelse(Z, 1, e/(1-e))) %>%   # for ATT
    summarize(ATThat = 
                sum(Z*wt*Y) / sum(Z*wt) -              # tx weighted mean
                sum((1-Z)*wt*Y) / sum((1-Z)*wt)) %>%   # co weighted mean
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
    summarize(ATThat = 
                sum(Z*wt*Y) / sum(Z*wt) -              # tx weighted mean
                sum((1-Z)*wt*Y) / sum((1-Z)*wt)) %>%   # co weighted mean
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

