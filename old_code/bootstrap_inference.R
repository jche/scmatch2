
# bootstrap methods

# # return single sample from dirichlet distribution with parameters alpha
# rdirichlet <- function(alpha) {
#   K <- length(alpha)
#   gammas <- rgamma(K, shape=alpha, scale=1)
#   return(gammas/sum(gammas))
# }


# bayesian bootstrap ------------------------------------------------------

#' Run a Bayesian bootstrap
#'
#' @param d dataset with weights for each tx/co unit
#' @param B number of bootstrap samples
#'
#' @return vector of bootstrapped ATT estimates
boot_bayesian <- function(d, B=100) {
  map_dbl(1:B,
          .progress = "Bootstrapping...",
          function(b) {
            d %>% 
              mutate(weights = weights * 
                       as.numeric(gtools::rdirichlet(1, alpha=rep(1,nrow(d))))) %>% 
                       # rdirichlet(alpha=rep(1,nrow(d)))) %>% 
              get_att_point_est()
          })
}


# bootstrap the covariates
boot_bayesian_covs <- function(d, covs, B=100) {
  map_dfr(1:B,
          .progress = "Bootstrapping...",
          function(b) {
            d %>% 
              mutate(weights = weights * 
                       as.numeric(gtools::rdirichlet(1, alpha=rep(1,nrow(d)))) *
                       nrow(d)) %>% 
              group_by(Z) %>% 
              summarize(across(all_of(covs),
                               ~sum(.*weights) / sum(weights))) %>% 
              summarize(across(all_of(covs),
                               ~last(.) - first(.)))
          })
}



#' Run a Bayesian bootstrap, corrected (?)
#'
#' @param d dataset with tx units and sc units
#' @param B number of bootstrap samples
#'
#' @return vector of bootstrapped ATT estimates
boot_bayesian2 <- function(d, B=100) {
  map_dbl(1:B,
          .progress = "Bootstrapping...",
          function(b) {
            set.seed(1)
            d %>% 
              mutate(
                boot_weights = rep(
                  gtools::rdirichlet(1, alpha=rep(1, nrow(d)/2)) %>% 
                    as.numeric(),
                  each=2),
                weights = weights * boot_weights) %>%
              get_att_point_est()
            
            d %>% 
              mutate(
                boot_weights = rep(
                  gtools::rdirichlet(1, alpha=rep(1, nrow(d)/2)) %>% 
                    as.numeric(),
                  each=2),
                tau_i = Z*Y - (1-Z)*Y) %>% 
              summarize(att = sum(boot_weights * tau_i)) %>% 
              pull(att)
          })
}


#' #' @param d dataset with tx units and sc units
#' #' @param B number of bootstrap samples
#' #'
#' #' @return vector of bootstrapped ATT estimates
#' boot_bayesian2_wild <- function(d, B=100) {
#'   map_dbl(1:B,
#'           .progress = "Bootstrapping...",
#'           function(b) {
#'             d %>% 
#'               mutate(
#'                 boot_weights = rep(
#'                   sample(c(-(sqrt(5)-1)/2, (sqrt(5)+1)/2), 
#'                          prob = c((sqrt(5)+1)/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5))),
#'                          size=nrow(d)/2, replace=T),
#'                   each=2),
#'                 weights = weights * boot_weights) %>%
#'               get_att_point_est()
#'           })
#' }


# directly compute T* = 1/sqrt(n1) sum(wt * (tau_i - tau))
boot_bayesian_redux <- function(d, B=100) {
  att_hat <- get_att_point_est(d)
  n <- nrow(d)
  n1 <- sum(d$Z)
  
  map_dbl(1:B,
          .progress = "Bootstrapping...",
          function(b) {
            d %>%
              mutate(weights_boot = as.numeric(gtools::rdirichlet(1, alpha=rep(1,nrow(d)))),
                     tau_i = Z*Y - (1-Z)*weights*Y) %>%
              summarize(att = sum(weights_boot * (tau_i - att_hat))) %>%
              pull(att)
          })
}


# input: df with weighted tx/co units
boot_bayesian_redux2 <- function(d, B=100) {
  map_dbl(1:B,
          .progress = "Bootstrapping...",
          function(b) {
            n <- nrow(d)
            n1 <- sum(d$Z)
            d %>% 
              mutate(weights_boot = as.numeric(gtools::rdirichlet(1, alpha=rep(1,nrow(d)))),
                     tau_i = Z*Y - (1-Z)*weights*Y) %>% 
              summarize(att = sum(weights_boot * (tau_i)) * n/n1) %>%
              pull(att)
          })
}


# input tx and aggregated co units
boot_bayesian_finalattempt <- function(d, B=100) {
  att_hat <- get_att_point_est(d)
  n1 <- sum(d$Z)
  n <- nrow(d)
  
  map_dbl(1:B,
          .progress = "Bootstrapping...",
          function(b) {
            d %>% 
              mutate(
                weights_boot = 
                  as.numeric(gtools::rdirichlet(1, alpha=rep(1,n))),
                  # sample(
                  #   c( -(sqrt(5)-1)/2, (sqrt(5)+1)/2 ),
                  #   prob = c( (sqrt(5)+1)/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5)) ),
                  #   replace = T, size = n),
                tau_i = Z*Y - (1-Z)*weights*Y) %>% 
              summarize(att = sum(weights_boot * (tau_i - att_hat)) / n1) %>%
              pull(att)
          })
}




# input tx and aggregated co units
#  - with wild bootstrap weights
boot_bayesian_finalattempt <- function(d, B=100) {
  att_hat <- get_att_point_est(d)
  n1 <- sum(d$Z)
  n <- nrow(d)
  
  map_dbl(1:B,
          .progress = "Bootstrapping...",
          function(b) {
            d %>% 
              mutate(
                weights_boot = 
                  as.numeric(gtools::rdirichlet(1, alpha=rep(1,n))),
                # sample(
                #   c( -(sqrt(5)-1)/2, (sqrt(5)+1)/2 ),
                #   prob = c( (sqrt(5)+1)/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5)) ),
                #   replace = T, size = n),
                tau_i = Z*Y - (1-Z)*weights*Y) %>% 
              summarize(att = sum(weights_boot * (tau_i - att_hat)) / n1) %>%
              pull(att)
          })
}


weighted_boot_naive <- function(d, wt_fun, B=100) {
  att_hat <- d %>%
    summarize(att = sum((2*Z-1) * weights * Y) / sum(Z)) %>%
    pull(att)
  
  map_dbl(1:B,
          .progress = "Bootstrapping...",
          function(b) {
            d %>% 
              mutate(
                weights_boot = wt_fun(n()),
                tau_i = (2*Z-1) * weights * Y) %>% 
              summarize(
                att = sum(weights_boot * 
                            (tau_i - sum(Z)/n()*att_hat)) / sum(Z) ) %>%
              pull(att)
          })
}
weighted_boot_norm <- function(d, wt_fun, B=100) {
  att_hat <- d %>%
    summarize(att = sum((2*Z-1) * weights * Y) / sum(Z)) %>%
    pull(att)
  
  map_dbl(1:B,
          .progress = "Bootstrapping...",
          function(b) {
            d %>% 
              mutate(
                weights_boot = wt_fun(n()),
                tau_i = (2*Z-1) * weights * Y) %>% 
              summarize(
                # n1 = sum(weights_boot*Z),
                # n = sum(weights_boot),
                att = sum(weights_boot * 
                            (tau_i - sum(weights_boot*Z)/n()*att_hat)) / 
                  sum(weights_boot*Z) ) %>%
              pull(att)
          })
}



# # this inputs tx and weighted co units
# wild_bootstrap <- function(d, B=100) {
#   # original ATT estimate
#   n <- nrow(d)
#   n1 <- sum(d$Z)
#   att_hat <- d %>%
#     summarize(att = sum((2*Z-1) * weights * Y) / n1) %>%
#     pull(att)
#   
#   map_dbl(1:B,
#           .progress = "Bootstrapping...",
#           function(b) {
#             d %>% 
#               mutate(
#                 weights_boot = sample(
#                   c( -(sqrt(5)-1)/2, (sqrt(5)+1)/2 ),
#                   prob = c( (sqrt(5)+1)/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5)) ),
#                   replace = T, size = n()),
#                 tau_i = (2*Z-1) * weights * Y) %>% 
#               summarize(
#                 att = sum(weights_boot * 
#                             (tau_i - n1/n*att_hat)) / n1 ) %>%
#               pull(att)
#           })
# }

# this inputs tx and sc units
wild_bootstrap <- function(d, B=100) {
  att_hat <- d %>%
    summarize(att = sum((2*Z-1) * Y) / sum(Z)) %>%
    pull(att)
  
  map_dbl(
    1:B,
    function(x) {
      d %>% 
        group_by(subclass) %>% 
        summarize(Yhat = sum((2*Z-1) * Y)) %>% 
        mutate(
          weights_boot = sample(
            c( -(sqrt(5)-1)/2, (sqrt(5)+1)/2 ),
            prob = c( (sqrt(5)+1)/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5)) ),
            replace = T, size = n())) %>% 
        summarize(
          att = sum(weights_boot * (Yhat - att_hat)) / n()
        ) %>% 
        pull(att)
    }
  )
}


# naive bootstrap ---------------------------------------------------------

#' Run naive bootstrap
#'
#' @param d original dataset
#' @param B number of bootstrap samples
#'
#' @return vector of bootstrapped ATT estimates
boot_naive <- function(d,
                       caliper,
                       metric,
                       rad_method,
                       dist_scaling,
                       est_method = "scm",
                       knn = 25,
                       B=50) {
  map_dfr(1:B,
          .progress = "Bootstrapping...",
          function(b) {
            # run full procedure
            res <- d %>% 
              sample_n(n(), replace=T) %>% 
              arrange(desc(Z)) %>% 
              mutate(id = 1:n()) %>% 
              get_cal_matches(caliper = caliper,
                              metric = metric,
                              rad_method = rad_method,
                              dist_scaling = dist_scaling,
                              est_method = est_method,
                              return = "all",
                              knn = 5)
            
            # get FSATT
            feasible_units <- attr(res, "adacalipers") %>% 
              filter(adacal <= caliper) %>% 
              pull(id)
            feasible_subclasses <- res %>% 
              filter(id %in% feasible_units) %>% 
              pull(subclass)
            
            return(tibble(FSATT = res %>% 
                            filter(subclass %in% feasible_subclasses) %>% 
                            get_att_point_est(),
                          SATT = res %>% 
                            get_att_point_est()))
          })
}












