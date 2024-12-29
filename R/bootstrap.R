
#' #' Old residual bootstrap code
#' #'
#' #' @param resids (desc)
#' #' @param B  (desc)
#' #' @param boot_mtd  (desc)
#' #' @param seed_addition  (desc)
#' #'
#' #' @return (desc)
#' #' @export
#' #'
#' boot_by_resids <-
#'   function(resids,
#'            B,
#'            boot_mtd,
#'            seed_addition = 123){
#'     T_star <- numeric(B)
#'     for (b in 1:B){
#'       set.seed(123 + seed_addition + b*13)
#'       n1 <- length(resids)
#'       # The implemented W is W(in the paper) / sqrt(n)
#'       if (boot_mtd=="Bayesian"){
#'         W = gtools::rdirichlet(1, alpha=rep(1,n1))
#'       }else if (boot_mtd=="wild"){
#'         W = sample(
#'           c( -(sqrt(5)-1)/2, (sqrt(5)+1)/2 ),
#'           prob = c( (sqrt(5)+1)/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5)) ),
#'           replace = T, size = n1) / n1
#'       }else if (boot_mtd == "sign"){
#'         W = sample(
#'           c(-1, 1),
#'           prob = c( 1/2, 1/2 ),
#'           replace = T, size = n1)
#'       }else if (boot_mtd == "naive-resid"){
#'         W_mat = rmultinom(1,n1, rep(1/n1,n1)) / 10
#'         W = c(W_mat)
#'       }
#'       T_star[b] = sum(resids * W)
#'     }
#'     return(T_star)
#'   }

#' Create a bootstrap function factory for residual resampling
#'
#' @description
#' Creates a specialized bootstrap function for resampling residuals using various methods.
#' The factory pattern allows creation of method-specific bootstrapping functions that can
#' be reused efficiently.
#'
#' @param boot_mtd The bootstrap method to use. Must be one of:
#'   \itemize{
#'     \item "Bayesian" - Uses Dirichlet distribution for weights
#'     \item "wild" - Uses wild bootstrap with Mammen weights
#'     \item "sign" - Uses random sign flipping
#'     \item "naive-resid" - Uses naive residual resampling
#'   }
#'
#' @return A function that performs the specified bootstrap resampling with parameters:
#'   \itemize{
#'     \item resids - Numeric vector of residuals to bootstrap
#'     \item B - Number of bootstrap iterations
#'     \item seed_addition - Value added to base seed for reproducibility (default: 123)
#'   }
#'
#' @examples
#' # Create bootstrap functions for different methods
#' bayesian_boot <- make_bootstrap("Bayesian")
#' wild_boot <- make_bootstrap("wild")
#'
#' # Generate some example residuals
#' resids <- rnorm(100)
#'
#' # Perform bootstrap with each method
#' bayesian_results <- bayesian_boot(resids, B = 1000)
#' wild_results <- wild_boot(resids, B = 1000)
#'
#' @export
make_bootstrap <- function(boot_mtd) {
  force(boot_mtd) # Ensure evaluation

  # Define the sampling function based on method
  sampler <- switch(boot_mtd,
                    "Bayesian" = function(n) {
                      gtools::rdirichlet(1, alpha = rep(1, n))
                    },
                    "wild" = function(n) {
                      sample(
                        c(-(sqrt(5)-1)/2, (sqrt(5)+1)/2),
                        prob = c((sqrt(5)+1)/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5))),
                        replace = TRUE,
                        size = n
                      ) / n
                    },
                    "sign" = function(n) {
                      sample(
                        c(-1, 1),
                        prob = c(1/2, 1/2),
                        replace = TRUE,
                        size = n
                      )
                    },
                    "naive-resid" = function(n) {
                      W_mat <- rmultinom(1, n, rep(1/n, n)) / 10
                      c(W_mat)
                    },
                    stop("Unknown bootstrap method:", boot_mtd)
  )

  #' Bootstrap Resampling Function
  #'
  #' @param resids Numeric vector of residuals to bootstrap
  #' @param B Number of bootstrap iterations
  #' @param seed_addition Value added to base seed for reproducibility
  #'
  #' @return Numeric vector of length B containing bootstrap statistics
  function(resids, B, seed_addition = 123) {
    n1 <- length(resids)
    T_star <- numeric(B)

    for(b in 1:B) {
      set.seed(123 + seed_addition + b*13)
      W <- sampler(n1)
      T_star[b] <- sum(resids * W)
    }

    T_star
  }
}

#' Create a bootstrap confidence interval calculator factory
#'
#' @description
#' Creates a specialized function for computing confidence intervals and standard errors
#' from bootstrap samples using various resampling methods.
#'
#' @param boot_mtd The bootstrap method to use. Must be one of:
#'   \itemize{
#'     \item "Bayesian" - Uses Dirichlet distribution for weights
#'     \item "wild" - Uses wild bootstrap with Mammen weights
#'     \item "naive-resid" - Uses naive residual resampling
#'   }
#'
#' @return A function that computes confidence intervals and standard errors with parameters:
#'   \itemize{
#'     \item resids - Numeric vector of residuals to bootstrap
#'     \item mean_est - Point estimate (typically mean) around which to construct intervals
#'     \item B - Number of bootstrap iterations
#'     \item seed_addition - Value added to base seed for reproducibility (default: 123)
#'   }
#'
#' @return A list containing:
#'   \itemize{
#'     \item ci_lower - Lower bound of confidence interval
#'     \item ci_upper - Upper bound of confidence interval
#'     \item sd - Bootstrap standard error
#'   }
#'
#' @examples
#' # Create bootstrap CI calculator
#' bayesian_ci <- make_bootstrap_ci("Bayesian")
#'
#' # Generate example data
#' resids <- rnorm(100)
#' mean_est <- mean(resids)
#'
#' # Calculate CIs and SE
#' results <- bayesian_ci(resids, mean_est, B = 1000)
#' print(c(results$ci_lower, results$ci_upper, results$sd))
#'
#' @export
make_bootstrap_ci <- function(boot_mtd) {
  # First create the base bootstrap sampler
  bootstrap_sampler <- make_bootstrap(boot_mtd)

  #' Bootstrap Confidence Interval Calculator
  #'
  #' @param resids Numeric vector of residuals to bootstrap
  #' @param mean_est Point estimate around which to construct intervals
  #' @param B Number of bootstrap iterations
  #' @param seed_addition Value added to base seed for reproducibility
  #'
  #' @return List containing ci_lower, ci_upper, and sd
  function(resids, mean_est, B, seed_addition = 123) {
    # Validate input
    if (!is.numeric(resids) || !is.numeric(mean_est) || !is.numeric(B)) {
      stop("All inputs must be numeric")
    }

    # Generate bootstrap samples
    T_star <- bootstrap_sampler(resids, B, seed_addition)

    # Calculate confidence intervals and standard error
    list(
      ci_lower = mean_est - quantile(T_star, 0.975),
      ci_upper = mean_est - quantile(T_star, 0.025),
      sd = sd(T_star)
    )
  }
}


#' Main function: Estimate the variance from the bootstrap
#'
#' @param matches_table The data frame of the matched table
#' @param outcome Name of the outcome variable (default "Y")
#' @param treatment Name of the treatment variable (default "Z")
#' @param var_weight_type The way that cluster variances are averaged
#' "num_units": weight by number of units in the subclass
#' "ess_units": weight effective size of units in the subclass
#' "uniform: weight each cluster equally
#' @return A tibble with SE, sigma_hat, N_T, and N_C_tilde
#' @export
boot_SE <- function(
    matches_table,
    outcome = "Y",
    treatment = "Z",
    B = 100,
    boot_mtd = "sign"){

  if ( is.csm_matches(matches_table) ) {
    matches_table <- full_unit_table(matches_table)
  }

  # Step 1: give uniform weight
  # get_att_ests(matches_table, outcome = "Y")
  nrow(matches_table %>% filter(Z==0) %>% distinct(id))
  matches_table_weights <-
    matches_table %>%
    group_by(subclass) %>%
    mutate(n_c = n() - 1,  # Number of control units (Z == 0) in the subclass
           weights_unif = ifelse(Z == 1, 1, 1 / n_c)) %>%
    ungroup() %>%
    rename(weights_SCM = "weights")

  # Step 2: get residuals
  ATT_residuals <-
    matches_table_weights %>%
    group_by(subclass) %>%
    summarise(
      Y_Z1 = Y[Z == 1],  # Extract Y for the Z == 1 data
      uniform_weighted_sum_Y_Z0 = sum(Y[Z == 0] * weights_unif[Z == 0]),
      SCM_weighted_sum_Y_Z0 = sum(Y[Z == 0] * weights_SCM[Z == 0])
    ) %>%
    mutate(
      residual_unif = Y_Z1 - uniform_weighted_sum_Y_Z0,
      residual_SCM = Y_Z1 - SCM_weighted_sum_Y_Z0,
      )

  # Step 3: bootstrap residuals.
  # B = 100
  booted_unif_weight <- boot_by_resids(
    resids = ATT_residuals$residual_unif,
    B = B,
    boot_mtd = boot_mtd
  )

  booted_SCM_weight <- boot_by_resids(
    resids = ATT_residuals$residual_SCM,
    B = B,
    boot_mtd = boot_mtd
  )

  return(tibble(
    SE_unif_weight = sd(booted_unif_weight),
    SE_SCM_weight = sd(booted_SCM_weight)
  ))
}
