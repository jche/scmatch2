
#' Title
#'
#' @param resids
#' @param B
#' @param boot_mtd
#' @param seed_addition
#'
#' @return
#' @export
#'
#' @examples
boot_by_resids <-
  function(resids,
           B,
           boot_mtd,
           seed_addition = 123){
    T_star <- numeric(B)
    for (b in 1:B){
      set.seed(123 + seed_addition + b*13)
      n1 <- length(resids)
      # The implemented W is W(in the paper) / sqrt(n)
      if (boot_mtd=="Bayesian"){
        W = gtools::rdirichlet(1, alpha=rep(1,n1))
      }else if (boot_mtd=="wild"){
        W = sample(
          c( -(sqrt(5)-1)/2, (sqrt(5)+1)/2 ),
          prob = c( (sqrt(5)+1)/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5)) ),
          replace = T, size = n1) / n1
      }else if (boot_mtd == "sign"){
        W = sample(
          c(-1, 1),
          prob = c( 1/2, 1/2 ),
          replace = T, size = n1)
      }
      T_star[b] = sum(resids * W)
    }
    return(T_star)
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
