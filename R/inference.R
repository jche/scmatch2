
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
                       as.numeric(gtools::rdirichlet(1, alpha=rep(1,nrow(d)))) *
                       # rdirichlet(alpha=rep(1,nrow(d))) * 
                       nrow(d)) %>% 
              get_att_ests()
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
                       cal_method,
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
                              cal_method = cal_method,
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
                            get_att_ests(),
                          SATT = res %>% 
                            get_att_ests()))
          })
}












