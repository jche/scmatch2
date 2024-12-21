library(tidyverse)
library(mvtnorm)
source("R/distance.R")
source("R/sc.R")
source("R/matching.R")
source("R/estimate.R")
source("R/inference.R")
source("R/sim_data.R")
source("R/wrappers.R")
source("R/utils.R")
source("R/bootstrap.R")


get_se_AE <- function(preds_csm){
  # 1. Get debiased units; Get the subclasses
  preds_csm <- preds_csm %>%
    mutate(Y_bias_corrected = Y - hat_mu_0)
  # 2. Filter the controls
  #     and the subclasses with n_controls >= 2
  preds_csm_filtered <-
    preds_csm %>%
    filter(Z==F) %>%
    group_by(subclass) %>%
    filter(n() >= 2) %>%
    ungroup()
  # 3. For each filtered subclass,
  #   calculate the cluster residual s_j^2 = se(debiased_units)

  cluster_var_df <-
    preds_csm_filtered %>%
    group_by(subclass) %>%
    summarise(nj = n(),
              var_cluster = var(Y_bias_corrected))
  # Edit 10 Feb 2024: do not use
  # var_cluster = weighted_var(Y_bias_corrected, weights))

  # 4. Get the weighted average of s_j^2, weighted by n_j, number of units in the subclass
  weighted_var_df <- cluster_var_df %>%
    summarise(weighted_var = weighted.mean(var_cluster, w = nj))
  sigma_hat <- sqrt(weighted_var_df$weighted_var)

  # 5. calculate N_T and N_C
  N_T <- nrow(preds_csm %>% filter(Z==T))
  tmp <- preds_csm %>%
    filter(Z==F) %>%
    group_by(id) %>%
    summarise(w_i_sq = sum(weights^2))
  N_C_tilde <- N_T^2 / sum(tmp$w_i_sq)
  # 6. Calculate the variance of the estimator
  res <- sqrt((1/N_T + 1/N_C_tilde)) * sigma_hat
  return(res)
}

# Run
  deg_overlap = "low"
  if (deg_overlap == "low"){
    ctr_dist <- 0.5
  }else if (deg_overlap == "medium"){
    ctr_dist <- 0.3
  }else if (deg_overlap == "high"){
    ctr_dist <- 0.1
  }


{
  df_dgp <- gen_one_toy(ctr_dist=ctr_dist)

  # The V matrix
  dist_scaling <- df_dgp %>%
    summarize(across(starts_with("X"),
                     function(x) {
                       if (is.numeric(x)) 6 / (max(x) - min(x))
                       else 1000
                     }))

  matches_and_debiased_residuals<-
    get_matches_and_debiased_residuals(
      dgp_name="toy", df_dgp,
      dist_scaling, mu_model="linear",n_split=1)
  list2env(matches_and_debiased_residuals,
           envir = environment())

  se_AE <- get_se_AE(preds_csm)
  print(paste0("Avi-Eli standard error is: ",
               signif(se_AE,3) ))
}


set.seed(123)
R <- 100
N_Ts <- N_Cs <- numeric(R * 3)
ctr_dists <- c(0.5, 0.3, 0.1)
full_ctr_dists <- rep(c(0.5, 0.3, 0.1),each=R)
deg_overlaps <- rep(c("low","mid","high"),each=R)
for (ctr_dist in ctr_dists){
  for (i in 1:R) {
    print(i)
    df_dgp <- gen_one_toy(ctr_dist=ctr_dist)

    # The V matrix
    dist_scaling <- df_dgp %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) 6 / (max(x) - min(x))
                         else 1000
                       }))

    matches_and_debiased_residuals<-
      get_matches_and_debiased_residuals(
        dgp_name="toy", df_dgp,
        dist_scaling, mu_model="linear",n_split=1)
    list2env(matches_and_debiased_residuals,
             envir = environment())

    list2env(calc_N_T_N_C(preds_csm),
             envir = environment())
    N_Ts[i] <- N_T
    N_Cs[i] <- N_C_tilde
    se_AE <- get_se_AE(preds_csm)
    # print(paste0("Avi-Eli standard error is: ",
    #              signif(se_AE,3) ))
  }
}

df_N_T_N_C <-
  data.frame(N_Ts=N_Ts,N_Cs=N_Cs,deg_overlaps=deg_overlaps)

FNAME <-paste0("data/outputs/sim_toy_results/","A_E_N_C_N_T.csv")
write_csv(df_N_T_N_C, FNAME, append=T)
