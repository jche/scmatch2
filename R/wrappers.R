# R/wrappers.R
# wrapper functions of other packages: an adaptive implementation of
# wide-use basic functions for convenience
#
# This includes a set of functions that simply estimates and returns
# the att as a single number, useful for the simulation studies.



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





get_att_diff <- function(d) {
  if ( is.csm_matches( d ) ) {
    d <- result_table(d, "sc" )
  }

  d %>%
    group_by(Z) %>%
    summarize(mn = mean(Y) ) %>%
    mutate( Z = as.numeric( Z != 0 ) ) %>%
    pivot_wider(names_from=Z, names_prefix="Y", values_from=mn) %>%
    mutate(ATT = Y1 - Y0 ) %>%
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
    filter(Z==1) %>%
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
                        scaling,
                        rad_method = "adaptive",
                        est_method = "scm") {
  preds_csm <- get_cal_matches(
    df = d,
    metric = metric,
    scaling = scaling,
    rad_method = rad_method,
    est_method = est_method )


  if (F) {
    # when return = "all":
    #  - in hain simulation, see that you get some
    #    very extreme estimates from units with 1 or 2 matched controls...
    result_table( preds_csm, "sc" ) %>%
      group_by(subclass) %>%
      summarize(n = n(),
                est = sum(weights*Y*Z) - sum(weights*Y*(1-Z))) %>%
      ggplot(aes(x=n-1, y=est)) +
      geom_point()
  }

  # output average number of matches per unit
  if (F) {
    n_matches <- preds_csm$matches %>%
      map_dbl(~nrow(.)-1)
    print(paste("Average number of co units per tx unit:", mean(n_matches)))
    p <- ggplot(tibble(x=n_matches)) +
      geom_histogram(aes(x), color="black")
    print(p)
  }

  get_att_diff( preds_csm )

}



#' Perform Matching Based on Specified Type
#'
#' This function performs matching on a dataset based on the specified `matching_type`.
#' It supports different matching methods and configurations, such as fixed-radius or
#' k-nearest-neighbor methods.
#'
#' @param matching_type A string specifying the type of matching. Supported types are:
#'   - `"maximum_fixed_scm"`: Uses the maximum metric with fixed radius and SCM estimation.
#'   - `"euclidean_knn"`: Uses the Euclidean metric with k-nearest neighbors.
#' @param df_dgp A data frame containing the dataset to be matched. The dataset should include
#'   treatment and covariate columns.
#' @param scaling A numeric value or vector used for scaling covariates during the matching process.
#'   Defaults to `1`.
#' @param k A numeric value to specify the number of nearest neighbors in knn matching
#'   Defaults to `8`
#'
#' @return A list containing the matched dataset and associated metrics. The specific structure of
#' the result depends on the `matching_type` and the underlying matching method used.
#'
#' @details
#' The function abstracts the matching process, allowing for different types of matching algorithms
#' to be applied to the same dataset. The input `matching_type` determines the matching method,
#' metric, and additional parameters.
#'
#' @examples
#' # Example for 'maximum_fixed_scm' matching type
#' test_df <- data.frame(
#'   Z = c(1, 0, 0, 0, 1),
#'   X = c(0, 0.5, 0.8, 3, 1.6)
#' )
#' scaling <- 1
#' result <- get_matches(
#'   matching_type = "maximum_fixed_scm",
#'   df_dgp = test_df,
#'   scaling = scaling
#' )
#'
#' # Example for 'euclidean_knn' matching type
#' test_df <- data.frame(
#'   Z = c(1, 0, 0, 0, 1),
#'   X1 = c(0, 0.5, 0.8, 3, 1.6),
#'   X2 = c(0, 0, 0, 0, 0)
#' )
#' result <- get_matches(
#'   matching_type = "euclidean_knn",
#'   df_dgp = test_df,
#'   scaling = scaling
#' )
#'
#' @export
get_matches <- function(matching_type,
                        df_dgp,
                        scaling,
                        k = 8) {
  if (matching_type == "maximum_fixed_scm") {
    df_dgp_with_matches <- get_cal_matches(
      df = df_dgp,
      metric = "maximum",
      scaling = scaling,
      rad_method = "fixed",
      est_method = "scm"
    )

  } else if (matching_type == "euclidean_knn") {
    df_dgp_with_matches <- get_cal_matches(
      df = df_dgp,
      covs = starts_with("X"),
      treatment = "Z",
      scaling = 1,
      metric = "euclidean",
      rad_method = "knn",
      k = k
    )
  } else {
    stop("Invalid matching_type. Must be 'maximum_fixed_scm' or 'euclidean_knn'.")
  }
  df_dgp_with_matches <- full_unit_table(
    df_dgp_with_matches,
    nonzero_weight_only = FALSE
  )
  return(df_dgp_with_matches)
}



# Code to implement CEM matching via the MatchIt package
# Allows SCM within cells if desired.
get_cem_matches <- function(
    df,
    covs = get_x_vars(df),
    Z_FORMULA = as.formula(paste0("Z~",
                                  paste0(grep("^X", names(df), value=T),
                                         collapse="+"))),
    num_bins,
    est_method = c("average", "scm"),
    return = c("sc_units", "agg_co_units", "all"), warn = TRUE ) {

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

  # calculate scaling
  scaling = df %>%
    summarize(across(all_of(covs),
                     function(x) {
                       if (is.numeric(x)) num_bins / (max(x) - min(x))
                       else 1000
                     }))

  # estimate outcomes within cells
  scweights <- est_weights(matched_gps = cem_matched_gps,
                           covs=covs,
                           scaling = scaling,
                           est_method = est_method,
                           metric = "maximum")

  # aggregate outcomes as desired
  m.data <- switch(return,
                   sc_units     = agg_sc_units(scweights),
                   agg_co_units = agg_co_units(scweights),
                   all          = bind_rows(scweights))

  unmatched_units <- setdiff(df %>% filter(Z==1) %>% pull(id),
                             m.data %>% filter(Z==1) %>% pull(id))
  if (warn && length(unmatched_units) > 0) {
    warning(glue::glue("Dropped the following treated units from data:
                        \t {paste(unmatched_units, collapse=\", \")}"))
  }

  attr(m.data, "feasible_units") <- m.data %>%
    filter(Z==1) %>%
    pull(id)

  return(m.data)
}


# estimates ATT!
get_att_cem <- function(d,
                        num_bins,
                        estimand = c("ATT", "CEM-ATT"),
                        est_method = "average",
                        extrapolate = TRUE ) {
  estimand <- match.arg(estimand)

  preds_feasible <- get_cem_matches(
    df = d,
    num_bins = num_bins,
    est_method = est_method,
    warn=FALSE )
  # print("Printing preds_feasible")
  # print(preds_feasible)
  # get ATT estimate:
  att_feasible <- get_att_diff( preds_feasible )

  if (estimand == "CEM-ATT") {
    return(att_feasible)
  }

  # Extrapolate to units not in bins?
  if (extrapolate && ( length(attr(preds_feasible, "feasible_units")) < sum(d$Z))) {
    preds_infeasible <- d %>%
      filter(!Z | !(id %in% attr(preds_feasible, "feasible_units"))) %>%
      get_cal_matches(.,
                      metric = "maximum",
                      rad_method = "1nn",
                      est_method = est_method,
                      scaling = d %>%
                        summarize(across(starts_with("X"),
                                         function(x) {
                                           if (is.numeric(x)) num_bins / (max(x) - min(x))
                                           else 1000
                                         })) )
    att_infeasible <- get_att_diff( preds_infeasible )

    return((att_feasible * sum(preds_feasible$Z) +
              att_infeasible * sum(preds_infeasible$Z)) / sum(d$Z))
  }
  return(att_feasible)
}




get_att_1nn <- function(d, scaling) {

  preds_1nn <- get_cal_matches(
    df = d,
    metric = "maximum",
    rad_method = "1nn",
    est_method = "average",
    scaling = scaling )

  # get ATT estimate:
  get_att_diff( preds_1nn )
}





run_all_methods <- function( df, skip_slow = TRUE, num_bins = 5, extrapolate = TRUE ) {

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

  one_nn <- get_att_1nn( df, scaling = 1/num_bins )

  if ( skip_slow ) {
    tmle = NA
    aipw = NA
  } else {
    # TODO: What are these libraries?  They don't seem to be defined.
    SL.library1 = c("SL.glm", "tmle.SL.dbarts2", "SL.glmnet")

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

  csm <- CSM:::get_att_csm(df, scaling = 1/num_bins, est_method="scm")

  cem <- CSM:::get_att_cem(df, num_bins=num_bins, est_method="scm",
                           extrapolate = extrapolate )

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

