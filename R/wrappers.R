
# wrapper functions:
#  an adaptive implementation of wide-use basic
#  functions for convenience

require(optweight)
require(dbarts)
require(tmle)
require(AIPW)
require(mvtnorm)


#' SuperLearner fitting
#'
#' @param df_to_fit The data frame for whom we fit the SuperLearner.
#' Should contain both X and Y
#' @param X_names A vector of strings, specifying column names of covariates
#' @param Y_name A string variable specifying
#' @param SL.library A string to tell SuperLearner which learner to use
#'
#' @return A SuperLearner object
#' @export
#'
#' @examples get_SL_fit <- function(df_to_fit=df_to_fit_i,X_names=c("X1","X2"),Y_name="Y", SL.library="SL.lm")
get_SL_fit <- function(df_to_fit,
                       X_names,
                       Y_name,
                       SL.library){
  SuperLearner(Y = df_to_fit[,Y_name,drop=T],
               X = df_to_fit[,X_names],
               SL.library = SL.library)
}

#' Get prediction from the fitted SuperLearner object
#'
#' @param SL_fit A fitted SuperLearner object
#' @param df_pred The testing dataset
#' @param X_names A vector of strings, specifying column names of covariates
#'
#' @return A matrix of predictions. Column is the SL.lib. Rows are N_data
#' @export
#'
#' @examples mock_SL_fit <- create_mock_SL_fit()
#' mock_df_test <- create_mock_df_test()
#' X_names <- c("X1","X2")
#' result <- get_SL_pred(mock_SL_fit,
#'                       mock_df_test,
#'                       X_names)
get_SL_pred <-
  function(SL_fit, df_test, X_names){
    SL_pred_obj <- predict(SL_fit,
                           newdata = df_test[,X_names])
    return(SL_pred_obj$library.predict)
  }


gen_one_toy <- function(ctr_dist = 0.5){
  gen_df_adv(
    nc=500,
    nt=100,
    f0_sd = 0.5,
    tx_effect_fun = function(X1, X2) {3*X1+3*X2},
    f0_fun = function(x,y) {
      matrix(c(x,y), ncol=2) %>%
        dmvnorm(mean = c(0.5,0.5),
                sigma = matrix(c(1,0.8,0.8,1), nrow=2)) * 20   # multiply for more slope!
    },
    ctr_dist = ctr_dist)
}


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
  trt_var_name <- all.vars(form)[1]
  trt_var <- d[[trt_var_name]]

  if (!all(trt_var %in% c(F,T) ) ){
    stop("Treatment variable must be either 0, 1 or FALSE, TRUE")
  }

  ## Fit model
  m_bal <- optweight(form,
                     data = d,
                     tols = tols,
                     estimand = "ATT")

  ### output ATT estimate
  d <- d %>%
    ungroup() %>%
    mutate(wt = m_bal$weights) %>%
    group_by(Z) %>% # Change this step
    summarize(Y_wtd = weighted.mean(Y,wt))

  att <- diff(d$Y_wtd)
  return(att)
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
                   select({{covs}}),
                 verbose=F)

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
                   select({{covs}}),
                 verbose=F)

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
                        metric = "maximum",
                        dist_scaling,
                        rad_method = "adaptive",
                        est_method = "scm") {
  preds_csm <- get_cal_matches(
    df = d,
    metric = metric,
    dist_scaling = dist_scaling,
    rad_method = rad_method,
    est_method = est_method,
    return = "sc_units",
    knn = 25)

  if (F) {
    # when return = "all":
    #  - in hain simulation, see that you get some
    #    very extreme estimates from units with 1 or 2 matched controls...
    preds_csm %>%
      group_by(subclass) %>%
      summarize(n = n(),
                est = sum(weights*Y*Z) - sum(weights*Y*(1-Z))) %>%
      ggplot(aes(x=n-1, y=est)) +
      geom_point()
  }

  # output average number of matches per unit
  if (F) {
    n_matches <- attr(preds_csm, "scweights") %>%
      map_dbl(~nrow(.)-1)
    print(paste("Average number of co units per tx unit:", mean(n_matches)))
    p <- ggplot(tibble(x=n_matches)) +
      geom_histogram(aes(x), color="black")
    print(p)
  }

  # get ATT estimate:
  preds_csm %>%
    group_by(Z) %>%
    summarize(Y = mean(Y)) %>%   # return = "sc_units" so this is okay
    pivot_wider(names_from=Z, names_prefix="Y", values_from=Y) %>%
    mutate(ATT = YTRUE - YFALSE) %>%
    pull(ATT)
}

# estimates ATT!
get_att_cem <- function(d,
                        num_bins,
                        estimand = c("ATT", "CEM-ATT"),
                        est_method = "average") {
  estimand <- match.arg(estimand)

  preds_feasible <- get_cem_matches(
    df = d,
    num_bins = num_bins,
    est_method = est_method,
    return = "sc_units")

  # get ATT estimate:
  att_feasible <- preds_feasible %>%
    group_by(Z) %>%
    summarize(Y = mean(Y)) %>%   # return = "sc_units" so this is okay
    pivot_wider(names_from=Z, names_prefix="Y", values_from=Y) %>%
    mutate(ATT = YTRUE - YFALSE) %>%
    pull(ATT)

  if (estimand == "CEM-ATT") {
    return(att_feasible)
  }

  if (length(attr(preds_feasible, "feasible_units")) < sum(d$Z)) {
    preds_infeasible <- d %>%
      filter(!Z | !(id %in% attr(preds_feasible, "feasible_units"))) %>%
      get_cal_matches(.,
                      metric = "maximum",
                      rad_method = "1nn",
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


get_att_1nn <- function(d, dist_scaling) {
  preds_1nn <- get_cal_matches(
    df = d,
    metric = "maximum",
    rad_method = "1nn",
    est_method = "average",
    dist_scaling = dist_scaling,
    return = "sc_units",
    cem_tx_units = d %>%
      filter(Z==1) %>%
      pull(id)
  )

  # get ATT estimate:
  preds_1nn %>%
    group_by(Z) %>%
    summarize(Y = mean(Y)) %>%   # return = "sc_units" so this is okay
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

