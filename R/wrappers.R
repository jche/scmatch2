
# wrapper functions of other packages: an adaptive implementation of
# wide-use basic functions for convenience
#
# This includes a set of functions that simply estimates and returns
# the att as a single number, useful for the simulaton studies.



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
  require(tmle)
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
get_SL_pred <- function(SL_fit, df_test, X_names){
  require(tmle)

  SL_pred_obj <- predict(SL_fit,
                         newdata = df_test[,X_names])
  return(SL_pred_obj$library.predict)
}

#' Caliper Synthetic Matching
#'
#' This function implements (adaptive) radius matching with optional
#' synthetic step on the resulting sets of controls.
#'
#' @param df The data frame to be matched
#' @param metric A string specifying the distance metric
#' @param caliper A numeric specifying the caliper size
#' @param rad_method adaptive caliper, fixed caliper, only 1nn caliper
#' @param est_method A string specifying the estimation method
#' @param return A string specifying what to return
#' @param dist_scaling A vector of scaling constants for covariates
#' @param ... Additional arguments
#'
#' @return df with a bunch of attributes.
#' @export
get_cal_matches <- function( df,
                             covs = starts_with("X"),
                             treatment = "Z",
                             metric = c("maximum", "euclidean", "manhattan"),
                             caliper = 1,
                             rad_method = c("adaptive", "fixed", "1nn"),
                             est_method = c("scm", "scm_extrap", "average"),
                             return = c("sc_units", "agg_co_units", "all"),
                             dist_scaling = df %>%
                               summarize(across(all_of(covs),
                                                function(x) {
                                                  if (is.numeric(x)) 1/sd(x)
                                                  else 1000
                                                })),
                             ...) {
  metric <- match.arg(metric)
  rad_method <- match.arg(rad_method)
  est_method <- match.arg(est_method)
  return <- match.arg(return)
  args <- list(...)
  df$id <- 1:nrow(df)

  ### use rad_method: generate matches
  # get caliper matches
  scmatches <- df %>%
    gen_matches(
      covs = covs,
      treatment = treatment,
      scaling = dist_scaling,
      metric = metric,
      caliper = caliper,
      rad_method = rad_method,
      ...)

  # scmatches$matches is a list of length ntx.
  #   Each element is a data frame of matched controls

  ### use est_method: scm or average
  scweights <- est_weights(df,
                           covs=covs,
                           matched_gps = scmatches$matches,
                           dist_scaling = dist_scaling,
                           est_method = est_method,
                           metric = metric)

  ### use return: sc units, aggregated weights per control, all weights
  m.data <- switch(return,
                   sc_units     = agg_sc_units(scweights),
                   agg_co_units = agg_co_units(scweights),
                   all          = bind_rows(scweights))

  unmatched_units <- setdiff(df %>% filter(.data[[treatment]]==1) %>% pull(id),
                             m.data %>% filter(.data[[treatment]]==1) %>% pull(id))
  if (length(unmatched_units) > 0) {
    if ( length( unmatched_units <= 20 ) ) {
      warning(glue::glue("Dropped the following treated units from data:
                        \t {paste(unmatched_units, collapse=\", \")}"))
    } else {
      warning(glue::glue("Dropped {length(unmatched_units) treated units from data.") )
    }
  }

  # aggregate results
  res <- m.data

  # store information about scaling and adaptive calipers
  adacalipers_df <- tibble(
    id = df %>% filter(Z == 1) %>% pull(id) %>% as.character(),
    adacal = scmatches$adacalipers)
  attr(res, "scaling") <- dist_scaling
  attr(res, "adacalipers") <- adacalipers_df

  # store information about feasible units/subclasses
  feasible_units <- adacalipers_df %>%
    filter(adacal <= caliper) %>%
    pull(id)
  attr(res, "unmatched_units") <- unmatched_units
  attr(res, "feasible_units")  <- feasible_units
  attr(res, "feasible_subclasses") <- res %>%
    filter(id %in% feasible_units) %>%
    pull(subclass)

  # keep a bunch of data around, just in case
  attr(res, "scweights") <- scweights
  attr(res, "dm") <- scmatches$dm
  attr(res, "dm_uncapped") <- scmatches$dm_uncapped

  return(res)
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
  m_bal <- optweight::optweight(form,
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
  require( dbarts )

  m_bart <- bart(x.train = d %>%
                   filter(Z==0) %>%
                   select( all_of(covs) ),
                 y.train = d %>%
                   filter(Z==0) %>%
                   pull(Y),
                 x.test = d %>%
                   filter(Z==1) %>%
                   select(all_of(covs)),
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

  require( dbarts )

  m_bart <- bart(x.train = d %>%
                   select(all_of(covs)),
                 y.train = d$Z %>% as.numeric(),
                 x.test = d %>%
                   select(all_of(covs)),
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

  require( tmle )

  tmle <- tmle(Y = d$Y,
               A = as.numeric(d$Z),
               W = d %>%
                 select(all_of(covs)),
               Q.SL.library = Q.SL.library,
               g.SL.library = g.SL.library,
               V.Q = 5, V.g = 5, V.Delta = 5, V.Z = 5)

  tmle$estimates$ATT$psi
}

get_att_aipw <- function(d,
                         covs,
                         Q.SL.library,
                         g.SL.library) {

  require( AIPW )
  aipw <- AIPW$
    new(Y = d$Y,
        A = d$Z,
        W = d %>%
          select(all_of(covs)),
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



# Code to implement CEM matching via the MatchIt package
get_cem_matches <- function(
    df,
    covs = get_x_vars(df),
    Z_FORMULA = as.formula(paste0("Z~",
                                  paste0(grep("^X", names(df), value=T),
                                         collapse="+"))),
    num_bins,
    est_method = c("average", "scm"),
    return = c("sc_units", "agg_co_units", "all")) {

  est_method <- match.arg(est_method)
  return <- match.arg(return)
  # print(Z_FORMULA)
  # print(head(df))
  m.out3 <- MatchIt::matchit(
    Z_FORMULA,
    data = df,
    method = "cem",
    estimand = "ATT",
    grouping = NULL,   # exact match factor covariates
    cutpoints = num_bins,
    k2k = F)   # not 1:1 matching
  m.data3 <- MatchIt::match.data(m.out3)
  # print(glue("Number of tx units kept: {sum(m.data3$Z)}"))

  # transform m.data3 into our format
  co_units <- m.data3 %>%
    filter(Z == 0)
  tx_units <- m.data3 %>%
    filter(Z == 1) %>%
    mutate(subclass_new = 1:n())  # give each tx unit its own subclass,
  # instead of letting multiple tx units share subclass

  # generate list of dfs, tx unit in first row
  # TODO: this is terribly slow, for no good reason...
  cem_matched_gps <- tx_units %>%
    split(f = ~subclass_new) %>%
    map(function(d) {
      d %>%
        bind_rows(
          co_units %>%
            filter(subclass == d$subclass) %>%
            mutate(subclass_new = d$subclass_new))
    }) %>%
    map(~.x %>%
          mutate(subclass = subclass_new) %>%
          select(-subclass_new, -weights))

  # estimate outcomes within cells
  scweights <- est_weights(df,
                           covs=names(df),
                           matched_gps = cem_matched_gps,
                           dist_scaling = df %>%
                             summarize(across(all_of(covs),
                                              function(x) {
                                                if (is.numeric(x)) num_bins / (max(x) - min(x))
                                                else 1000
                                              })),
                           est_method = est_method,
                           metric = "maximum")

  # aggregate outcomes as desired
  m.data <- switch(return,
                   sc_units     = agg_sc_units(scweights),
                   agg_co_units = agg_co_units(scweights),
                   all          = bind_rows(scweights))

  unmatched_units <- setdiff(df %>% filter(Z==1) %>% pull(id),
                             m.data %>% filter(Z==1) %>% pull(id))
  if (length(unmatched_units) > 0) {
    warning(glue::glue("Dropped the following treated units from data:
                        \t {paste(unmatched_units, collapse=\", \")}"))
  }

  attr(m.data, "feasible_units") <- m.data %>%
    filter(Z) %>%
    pull(id)

  return(m.data)
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
  # print("Printing preds_feasible")
  # print(preds_feasible)
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





run_all_methods <- function( df, skip_slow = TRUE ) {

  bal <- CSM:::get_att_bal(df,
                           form=as.formula('Z ~ X1+X2'),
                           tols=c(0.01, 0.01))

  bal_int <- CSM:::get_att_bal(df,
                               form=as.formula('Z ~ X1*X2'),
                               tols=c(0.01, 0.01, 0.1))

  or_lm <- CSM:::get_att_or_lm(df, form=as.formula('Y ~ X1*X2'))

  or_bart <- CSM:::get_att_or_bart(df, covs=c("X1", "X2"))

  lm_ps <- CSM:::get_att_ps_lm(df, form=as.formula('Z ~ X1*X2'))

  ps_bart <- CSM:::get_att_ps_bart(df, covs=c("X1", "X2"))

  one_nn <- get_att_1nn( df, dist_scaling = 1/5 )

  # TODO: What are these libraries?  They don't seem to be defined.
  SL.library1 = c("SL.glm", "tmle.SL.dbarts2", "SL.glmnet")

  if ( skip_slow ) {
    tmle = NA
    aipw = NA
  } else {
    tmle <- CSM:::get_att_tmle(df %>%
                                 mutate(X3 = X1*X2),
                               covs=c("X1","X2","X3"),
                               Q.SL.library = SL.library1,
                               g.SL.library = SL.library1)

    aipw <- CSM:::get_att_aipw(df %>%
                                 mutate(X3 = X1*X2),
                               covs=c("X1","X2","X3"),
                               Q.SL.library = SL.library1,
                               g.SL.library = SL.library1)
  }
  csm <- CSM:::get_att_csm(df, dist_scaling = 1/5, est_method="scm")

  cem <- CSM:::get_att_cem(df, num_bins=5, est_method="scm")

  tibble(method = c("bal", "bal_int", "or_lm", "or_bart",
                    "lm_ps", "ps_bart", "tmle", "aipw",
                    "csm", "cem", "1nn"),
         att = c(bal, bal_int, or_lm, or_bart,
                 lm_ps, ps_bart, tmle, aipw,
                 csm, cem, one_nn) )
}

# test these --------------------------------------------------------------

if (FALSE) {
  require(tidyverse)
  require( CSM )

  # df <- CSM:::gen_df_adv2d(nc=1000, nt=50,
  #                    tx_effect=0.2,
  #                    sd=0.03,
  #                    effect_fun = function(x,y) {
  #                      matrix(c(x,y), ncol=2) %>%
  #                        mvtnorm::dmvnorm(mean = c(0.5,0.5),
  #                                         sigma = matrix(c(1,0.8,0.8,1), nrow=2))
  #                    })

  # df <- CSM:::gen_df_full(nc=1000, nt=50, eps_sd=0,
  #                   tx_effect = function(x,y) {(0.5*(x+y))},
  #                   effect_fun = function(x,y) {x+y})

  df <- CSM:::gen_one_toy()

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


  run_all_methods( df )

}

