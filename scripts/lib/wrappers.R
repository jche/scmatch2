# scripts/wrappers.R
# wrapper functions of other packages: an adaptive implementation of
# wide-use basic functions for convenience
#
# This includes a set of functions that simply estimates and returns
# the att as a single number, useful for the simulation studies.
#
# Convention: functions that need covariates accept a formula of the form
#   Z ~ X1 + X2 + X3
# where the LHS is the treatment variable and the RHS lists the covariates.
# Use parse_form() to extract both from any such formula.


#' Extract treatment variable name and covariate names from a formula.
#'
#' @param form Formula of the form treatment ~ cov1 + cov2 + ...
#' @return Named list with elements:
#'   \item{treatment}{Character string: name of the treatment variable (LHS).}
#'   \item{covs}{Character vector: covariate names (RHS).}
#' @export
parse_form <- function(form) {
  vars <- all.vars(form)
  list(treatment = vars[1], covs = vars[-1])
}



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
    d <- result_table( d, "sc", outcome="Y" )
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
                            form) {
  require( dbarts )
  covs <- parse_form(form)$covs

  m_bart <- bart(x.train = d %>%
                   filter(Z==0) %>%
                   dplyr::select( all_of(covs) ),
                 y.train = d %>%
                   filter(Z==0) %>%
                   pull(Y),
                 x.test = d %>%
                   filter(Z==1) %>%
                   dplyr::select(all_of(covs)),
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
    mutate(e = CSM:::invlogit(predict(m_lm_ps, newdata=.)),
           wt = ifelse(Z, 1, e/(1-e))) %>%   # for ATT
    summarize(ATThat =
                sum(Z*wt*Y) / sum(Z*wt) -              # tx weighted mean
                sum((1-Z)*wt*Y) / sum((1-Z)*wt)) %>%   # co weighted mean
    pull(ATThat)
}

get_att_ps_bart <- function(d,
                            form) {

  require( dbarts )
  covs <- parse_form(form)$covs

  m_bart <- bart(x.train = d %>%
                   dplyr::select(all_of(covs)),
                 y.train = d$Z %>% as.numeric(),
                 x.test = d %>%
                   dplyr::select(all_of(covs)),
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
                         form,
                         Q.SL.library,
                         g.SL.library) {

  require( tmle )
  covs <- parse_form(form)$covs

  tmle <- tmle(Y = d$Y,
               A = as.numeric(d$Z),
               W = d %>%
                 dplyr::select(all_of(covs)),
               Q.SL.library = Q.SL.library,
               g.SL.library = g.SL.library,
               V.Q = 5, V.g = 5, V.Delta = 5, V.Z = 5)

  tmle$estimates$ATT$psi
}

get_att_aipw <- function(d,
                         form,
                         Q.SL.library,
                         g.SL.library) {

  require( AIPW )
  covs <- parse_form(form)$covs

  aipw <- AIPW$
    new(Y = d$Y,
        A = d$Z,
        W = d %>%
          dplyr::select(all_of(covs)),
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
                        est_method = "scm",
                        warn = FALSE ) {
  preds_csm <- get_cal_matches(
    data = d,
    treatment = "Z",
    metric = metric,
    scaling = scaling,
    rad_method = rad_method,
    est_method = est_method,
    warn = warn )


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
#' This function performs matching on a dataset based on the specified
#' `matching_type`. It supports different matching methods and
#' configurations, such as fixed-radius or k-nearest-neighbor methods.
#'
#' @param matching_type A string specifying the type of matching.
#'   Supported types are:
#'   - `"maximum_fixed_scm"`: Uses the maximum metric with fixed radius and SCM estimation.
#'   - `"euclidean_knn"`: Uses the Euclidean metric with k-nearest neighbors.
#' @param df_dgp A data frame containing the dataset to be matched.
#'   The dataset should include treatment and covariate columns.
#' @param scaling A numeric value or vector used for scaling
#'   covariates during the matching process. Defaults to `1`.
#' @param k A numeric value to specify the number of nearest neighbors
#'   in knn matching Defaults to `8`
#'
#' @return A list containing the matched dataset and associated
#'   metrics. The specific structure of the result depends on the
#'   `matching_type` and the underlying matching method used.
#'
#' @details The function abstracts the matching process, allowing for
#' different types of matching algorithms to be applied to the same
#' dataset. The input `matching_type` determines the matching method,
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
get_matches <- function( matching_type,
                         df_dgp,
                         scaling,
                         k = 8 ) {
  if (matching_type == "maximum_fixed_scm") {
    df_dgp_with_matches <- get_cal_matches(
      data = df_dgp,
      treatment="Z",
      metric = "maximum",
      scaling = scaling,
      rad_method = "fixed",
      est_method = "scm"
    )

  } else if (matching_type == "euclidean_knn") {
    df_dgp_with_matches <- get_cal_matches(
      data = df_dgp,
      covs = NULL,
      treatment = "Z",
      scaling = 1,
      metric = "euclidean",
      rad_method = "knn",
      k = k
    )
  } else {
    stop("Invalid matching_type. Must be 'maximum_fixed_scm' or 'euclidean_knn'.")
  }
  df_dgp_with_matches <- result_table(
    df_dgp_with_matches,
    nonzero_weight_only = FALSE
  )
  return(df_dgp_with_matches)
}



# Code to implement CEM matching via the MatchIt package
# Allows SCM within cells if desired.
get_cem_matches <- function(
    data,
    covs = CSM:::get_x_vars(data),
    Z_FORMULA = as.formula(paste0("Z~",
                                  paste0(grep("^X", names(data), value=T),
                                         collapse="+"))),
    num_bins,
    cutpoints = NULL,
    est_method = c("average", "scm"),
    return = c("sc_units", "agg_co_units", "all"), warn = TRUE ) {

  if ( is.character(Z_FORMULA)) {
    Z_FORMULA <- as.formula(paste0( "Z~",
                                    paste0(covs,collapse="+")) )
  }

  if ( is.null( cutpoints ) ) {
    cutpoints = num_bins
  }

  if ( !( "id" %in% names(data) ) ) {
    data <- data %>%
      mutate(id = as.character(1:n()))
  }


  est_method <- match.arg(est_method)
  return <- match.arg(return)
  # print(Z_FORMULA)
  # print(head(data))
  m.out3 <- MatchIt::matchit(
    Z_FORMULA,
    data = data,
    method = "cem",
    estimand = "ATT",
    grouping = NULL,   # exact match factor covariates
    cutpoints = cutpoints,
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
          dplyr::select(-subclass_new, -weights))

  # calculate scaling
  if ( is.list( cutpoints ) ) {
    scaling = map_dbl(cutpoints, ~ mean(diff(.x)))
  } else {
    scaling = data %>%
      summarize(across(all_of(covs),
                       function(x) {
                         if (is.numeric(x)) num_bins / (max(x) - min(x))
                         else 1000
                       }))
  }

  # estimate outcomes within cells
  treatment_var <- all.vars(Z_FORMULA)[1]
  scweights <- est_weights(matched_gps = cem_matched_gps,
                           covs=covs,
                           scaling = scaling,
                           est_method = est_method,
                           treatment = treatment_var,
                           metric = "maximum")

  # aggregate outcomes as desired
  m.data <- switch(return,
                   sc_units     = agg_sc_units(scweights,outcome="Y"),
                   agg_co_units = agg_co_units(scweights,outcome="Y"),
                   all          = bind_rows(scweights))

  unmatched_units <- setdiff(data %>% filter(Z==1) %>% pull(id),
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
    data = d,
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
                      treatment = "Z",
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
    data = d,
    treatment = "Z",
    metric = "maximum",
    rad_method = "1nn",
    est_method = "average",
    scaling = scaling )

  # get ATT estimate:
  get_att_diff( preds_1nn )
}


#' Estimate ATT using Causal Forests (grf package)
#'
#' @param d A data frame with treatment Z, outcome Y, and covariates
#' @param form Formula of the form Z ~ X1 + X2 + ... specifying treatment and covariates
#'
#' @return Estimated ATT
#' @export
get_att_causal_forest <- function(d, form) {
  if (!requireNamespace("grf", quietly = TRUE)) {
    stop("Package 'grf' is needed for this function. Please install it with: install.packages('grf')",
         call. = FALSE)
  }
  covs <- parse_form(form)$covs

  # Fit causal forest
  cf <- grf::causal_forest(
    X = as.matrix(d[, covs]),
    Y = d$Y,
    W = as.numeric(d$Z),
    num.trees = 2000
  )

  # Get average treatment effect on the treated
  # Predict CATE for all units
  tau_hat <- predict(cf)$predictions

  # Average over treated units only
  att <- mean(tau_hat[d$Z == 1])

  return(att)
}

#' Estimate ATT using TWANG (robust to logical/character covariates)
#'
#' @param d Data.frame with Y, Z, and covariates
#' @param form Propensity formula (e.g., Z ~ X1+X2+...)
#' @param stop.method TWANG stop method (default "es.mean")
#' @param n.trees GBM trees (default 5000)
#' @param interaction.depth GBM depth (default 3)
#' @param estimand "ATT" (default)
#' @return numeric(1) ATT or NA_real_ on failure
get_att_twang <- function(
    d, form,
    stop.method = "es.mean",
    n.trees = 5000L,
    interaction.depth = 3L,
    estimand = "ATT"
) {
  if (!requireNamespace("twang", quietly = TRUE)) {
    stop("Package 'twang' is needed. Install with install.packages('twang')", call. = FALSE)
  }

  # Ensure base type is data.frame (not tibble)
  d <- as.data.frame(d, stringsAsFactors = FALSE)

  # Ensure outcome is numeric
  if (!is.numeric(d$Y)) d$Y <- as.numeric(d$Y)

  # Ensure Z is 0/1 numeric
  if (is.logical(d$Z)) {
    d$Z <- as.integer(d$Z)
  } else if (is.factor(d$Z)) {
    d$Z <- as.integer(as.character(d$Z))
  } else {
    d$Z <- as.numeric(d$Z)
  }
  if (!all(d$Z %in% c(0, 1))) {
    warning("Treatment variable Z must be 0/1 after conversion")
    return(NA_real_)
  }
  if (length(unique(d$Z)) < 2) {
    warning("Only one treatment group present")
    return(NA_real_)
  }

  # Coerce all logical covariates to factor; characters to factor as well.
  is_logical <- vapply(d, is.logical, logical(1))
  if (any(is_logical)) {
    for (nm in names(d)[is_logical]) d[[nm]] <- factor(d[[nm]], levels = c(FALSE, TRUE))
  }
  is_char <- vapply(d, is.character, logical(1))
  if (any(is_char)) {
    for (nm in names(d)[is_char]) d[[nm]] <- factor(d[[nm]])
  }

  # Fit TWANG PS model
  ps_fit <- tryCatch(
    twang::ps(
      formula = form,
      data = d,
      estimand = estimand,
      stop.method = stop.method,
      n.trees = n.trees,
      interaction.depth = interaction.depth,
      verbose = FALSE
    ),
    error = function(e) {
      warning(sprintf("twang::ps failed: %s", e$message))
      return(NULL)
    }
  )
  if (is.null(ps_fit)) return(NA_real_)

  # Extract stabilized weights for chosen stop.method
  weights <- tryCatch(
    twang::get.weights(ps_fit, stop.method = stop.method),
    error = function(e) {
      warning(sprintf("get.weights failed: %s", e$message))
      return(NULL)
    }
  )
  if (is.null(weights) || !any(is.finite(weights))) {
    warning("Weights missing or non-finite")
    return(NA_real_)
  }

  # Weighted ATT: treated mean (weights ≈ 1) minus reweighted controls
  d$wt <- weights
  grp <- aggregate(list(Y = d$Y * d$wt, w = d$wt), by = list(Z = d$Z), FUN = sum)
  if (!all(c(0,1) %in% grp$Z) || any(grp$w == 0)) {
    warning("Invalid weights by group")
    return(NA_real_)
  }
  mu1 <- grp$Y[grp$Z == 1] / grp$w[grp$Z == 1]
  mu0 <- grp$Y[grp$Z == 0] / grp$w[grp$Z == 0]
  as.numeric(mu1 - mu0)
}

#' Estimate ATT using Kernel Balancing (kbal package)
#'
#' @param d A data frame with treatment Z, outcome Y, and covariates
#' @param form Formula of the form Z ~ X1 + X2 + ... specifying treatment and covariates
#' @param numdims Fixed number of dimensions to use (default: NULL for auto-selection)
#'
#' @return Estimated ATT
#' @export
get_att_kbal <- function(d, form, numdims = NULL) {
  if (!requireNamespace("kbal", quietly = TRUE)) {
    stop("Package 'kbal' is needed for this function. Please install it with: install.packages('kbal')",
         call. = FALSE)
  }
  covs <- parse_form(form)$covs

  # Check for valid Treatment variable (0/1 or T/F)
  # This ensures the test case for 'df_bad' passes
  if (!all(unique(d$Z) %in% c(0, 1, FALSE, TRUE))) {
    stop("Treatment variable must be either 0, 1 or FALSE, TRUE")
  }

  # Prepare data
  X <- as.matrix(d[, covs])
  Z <- as.numeric(d$Z)
  Y <- d$Y

  # Fit kernel balancing
  kbal_fit <- kbal::kbal(
    allx = X,
    treatment = Z,
    numdims = numdims,
    printprogress = FALSE  # Optional: suppress the long output
  )

  # kbal returns weights for ALL units (length = nrow(d))
  # For ATT: treated units have weight 1, control units get kbal weights
  weights <- ifelse(Z == 1, 1, kbal_fit$w)

  # Calculate weighted ATT
  d_weighted <- d %>%
    dplyr::mutate(wt = weights)

  att <- d_weighted %>%
    dplyr::group_by(Z) %>%
    dplyr::summarize(Y_wtd = weighted.mean(Y, wt)) %>%
    dplyr::summarize(att = diff(Y_wtd)) %>%
    dplyr::pull(att)

  return(att)
}
