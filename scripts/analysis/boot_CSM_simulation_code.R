
# Code to run the section5.4 simulation



#### Helper functions ####


get_df_scaling_from_dgp_name <- function(dgp_name,
                                          kang_true = F,
                                          toy_ctr_dist=0.5){
  if (dgp_name == "toy"){
    df_dgp <- gen_one_toy(ctr_dist=toy_ctr_dist)
    # Edit 4 Mar 2024:
    scaling <- 8
    # scaling <- df_dgp %>%
    #   summarize(across(starts_with("X"),
    #                    function(x) {
    #                      if (is.numeric(x)) 6 / (max(x) - min(x))
    #                      else 1000
    #                    }))
  }else if(dgp_name=="kang"){
    df_dgp <- gen_df_kang(n = 1000)
    if (kang_true){ # Replace X1 to X4 by V1 to V4
      df_dgp[, paste0("X",1:4)] <-
        df_dgp %>% select(starts_with("V"))
    }
    scaling <- df_dgp %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) 5 / (max(x) - min(x))
                         else 1000
                       }))
  }else if (dgp_name == "otsu"){
    df_dgp <- generate_one_otsu()
    scaling <- df_dgp %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) 6 / (max(x) - min(x))
                         else 1000
                       }))
  } else {
    stop("dgp_name must be toy or kang or otsu")
  }
  return(list(df_dgp=df_dgp,
              scaling=scaling))
}



get_matches_and_debiased_residuals <-
  function(dgp_name,
           df_dgp,
           scaling,
           mu_model,
           n_split=1){


    # TODO: What is the following for?  What is n_split for?  Not
    # following this code:

    preds_csm <- get_matches(dgp_name=dgp_name,
                             df_dgp=df_dgp,
                             scaling=scaling)

    #Assign group_label to df_dgp
    df_dgp_splitted <- assign_group_label(df_dgp,
                                          n_split = n_split) %>%
      mutate( id = as.character(id) )

    id_to_group_map <-
      df_dgp_splitted %>%
      select(id, group_label)

    # Append the same label to preds_csm through join
    preds_csm <-
      left_join(preds_csm$result,
                id_to_group_map %>% distinct(id, .keep_all=T),
                by="id")


    # Calculate pred_label
    pred_label_map <-
      get_pred_label_map(n=n_split)

    preds_csm <-
      left_join(preds_csm, pred_label_map, by="group_label")

    df_controls <-
      df_dgp_splitted %>% filter(Z==0)

    if (dgp_name == "toy" || dgp_name == "otsu") {
      X_names<- c("X1","X2")
      if (mu_model == "linear"){
        SL_lib <- "SL.lm"
      }else if (mu_model == "non-linear"){
        SL_lib <-  "SL.randomForest"
      }
    } else if (dgp_name == "kang") {
      if (mu_model == "kang_correct"){
        X_names<- c("V1","V2","V3","V4")
        SL_lib <- "SL.lm"
      }else if (mu_model == "linear"){
        X_names<- c("X1","X2","X3","X4")
        SL_lib <- "SL.lm"
      }
    }
    else {
      stop("dgp_name must be toy or kang")
    }


    # Model fitting using sample splitting
    # Nest data
    df_to_fit_nested <-
      df_controls %>%
      group_by(group_label) %>%
      nest() %>%
      arrange(group_label)
    models_n_split <-
      vector("list", length = n_split)
    for (i in 1:n_split){
      # i <- 1
      df_to_fit_i <- df_to_fit_nested$data[[i]]

      models_n_split[[i]]<-
        SL_fit <- get_SL_fit(df_to_fit=df_to_fit_i,
                             X_names = X_names,
                             Y_name = "Y",
                             SL.library = SL_lib)
      # get predictions:
      preds_csm[,paste0("hat_mu_0_pred_",i)] <-
        SL_pred  <- get_SL_pred(SL_fit=SL_fit,
                                df_test=preds_csm,
                                X_names=X_names)
    }

    # Finally, select the prediction
    preds_csm <- preds_csm %>%
      rowwise() %>%
      mutate(hat_mu_0 = get(paste0("hat_mu_0_pred_", pred_label)))


    ## Construct residuals
    # First, construct each \tilde \tau_i
    tmp0 <- preds_csm %>%
      mutate(Y_bias_corrected = Y - hat_mu_0) %>%
      group_by(subclass, Z) %>%
      summarize(mn = sum(Y_bias_corrected*weights))
    tmp <- tmp0 %>%
      group_by(subclass) %>%
      summarise(tilde_tau = last(mn)-first(mn), .groups="drop")
    tilde_tau = tmp$tilde_tau

    # Second, residual = \tilde \tau_i - mean(\tilde \tau_i)
    mean_tilde_tau <- mean(tilde_tau)
    tilde_tau_resids <- tilde_tau - mean_tilde_tau
    return(list(preds_csm=preds_csm,
                tilde_tau=tilde_tau,
                mean_tilde_tau=mean_tilde_tau,
                tilde_tau_resids=tilde_tau_resids))
  }



get_matches <- function(dgp_name, df_dgp, scaling){
  # TODO: WHY should the type of dgp affect how to get matches?  E.g.,
  # why is otsu different?  It should be a flag for kind of matching,
  # instead.

  if (dgp_name == "toy" || dgp_name == "kang") {
    df_dgp_with_matches <- get_cal_matches(
      df = df_dgp,
      metric = "maximum",
      scaling = scaling,
      rad_method = "fixed",
      est_method = "scm",
      return = "all"
    )
  } else if (dgp_name == "otsu") {
    df_dgp_with_matches <- get_NN_matches(df_dgp)

  } else {
    stop("dgp_name must be toy, kang, or otsu")
  }
  return(df_dgp_with_matches)
}



get_pred_label_map <- function(n){
  if (n == 1){
    return(tibble(group_label=1, pred_label=1))
  }else{
    return(tibble(group_label=1:n,
                  pred_label=c(2:n, 1) )) # index for cross-fitting
  }
}



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
  weighted_var <- function(x, wt) {
    n <- length(x)
    wm <- weighted.mean(x, wt)
    sum(wt * (x - wm)^2) * n / (n-1)
  }

  weighted_se <- function(x, wt) {
    sqrt(weighted_var(x, wt) / length(x))
  }
  cluster_var_df <-
    preds_csm_filtered %>%
    group_by(subclass) %>%
    summarise(nj = n(),
              var_cluster = var(Y_bias_corrected), .groups="drop")
  # var_cluster = weighted_var(Y_bias_corrected, weights))
  # 4. Get the weighted average of s_j^2, weighted by n_j, number of units in the subclass
  weighted_var_df <- cluster_var_df %>%
    summarise(weighted_var = weighted.mean(var_cluster, w = nj), .groups="drop")
  sigma_hat <- sqrt(weighted_var_df$weighted_var)

  # 5. calculate N_T and N_C
  list2env(calc_N_T_N_C(preds_csm),
           envir = environment())

  # 6. Calculate the variance of the estimator
  res <- sqrt((1/N_T + 1/N_C_tilde)) * sigma_hat
  return(res)
}




calc_N_T_N_C <- function(preds_csm){
  if ( is.csm_matches( preds_csm ) ) {
    preds_csm <- full_unit_table(preds_csm)
  }

  N_T <- nrow(preds_csm %>% filter(Z==T))
  tmp <- preds_csm %>%
    filter(Z==F) %>%
    group_by(id) %>%
    summarise(w_i = sum(weights), .groups="drop")
  N_C_tilde <- N_T^2 / sum(tmp$w_i^2)
  return(list(N_T = N_T,
              N_C_tilde = N_C_tilde ))
}


get_NN_matches <- function(df_dgp){
  library(MatchIt)

  m.out <- matchit(Z ~ X1 + X2,
                   data = df_dgp,
                   method = "nearest",
                   replace = TRUE,
                   mahvars = ~ X1 + X2,
                   ratio = 8)

  # # Obtain matched dataset
  # df_matched <- match.data(m.out)

  df_matched_ids <- data.frame(m.out$match.matrix)

  # # Inverting the treatment indicator
  # df_dgp$Z_inv <- ifelse(df_dgp$Z == 1, 0, 1)
  #
  # m.out2 <- matchit(Z_inv ~ X1 + X2,
  #                   data = df_dgp,
  #                   method = "nearest",
  #                   replace = TRUE,
  #                   mahvars = ~ X1 + X2,
  #                   ratio = 8)
  # df_matched_ids.2 <- data.frame(m.out2$match.matrix)
  #
  # df_matched_ids <- rbind(df_matched_ids,
  #                         df_matched_ids.2)

  id <- subclass <- weights <- c()

  for(i in 1:nrow(df_matched_ids)) {
    # Concatenate row name (treated id) with matched control ids
    # i <- 1
    row_ids <- c(as.integer(rownames(df_matched_ids)[i]),
                 as.integer(unlist(df_matched_ids[i, ])))
    id <- c(id, row_ids)
    # Add subclass
    M <- length(row_ids) - 1
    subclass <- c(subclass, rep(i, M + 1))
    weights <- c(weights, c(1, rep(1/M, M )) )
  }

  # Create the new dataframe
  new_df <- data.frame(id = id,
                       subclass = subclass,
                       weights)

  # Display the first few rows of the new dataframe
  head(new_df)

  preds_matching <-
    left_join(new_df, df_dgp%>%distinct(id, .keep_all=T), by = "id")
  return(preds_matching)
}

# Function: boot_bayesian
# input:
#   dgp_name: "toy" or "kang",
#   att0 = T or F: T means set true ATT to 0; F means get the true ATT from the
#   I = 100: number of MC runs
#   B = 250: number of bootstrap samples per MC run
#   mu_model="linear": string
#
# output: A tibble: id=1:I, name=dgp_name,
# boot_mtd = "Bayesian", att_true = att_true,
# att_est= att_est, att_est_debiased = att_debiased,
# lower=CI_lower,upper=CI_upper, covered=covered, sd_boot=sd_boot)
boot_bayesian <- function(dgp_name,
                          att0,
                          I=100,
                          B=250,
                          mu_model="linear",
                          n_split=1){
  boot_CSM(dgp_name,
           att0,
           I,
           B,
           mu_model,
           boot_mtd="Bayesian",
           n_split=n_split)
}

boot_wild <- function(dgp_name,
                      att0,
                      I=100,
                      B=250,
                      mu_model="linear",
                      n_split=1){
  boot_CSM(dgp_name,
           att0,
           I,
           B,
           mu_model,
           boot_mtd="wild",
           n_split=n_split)
}

boot_by_resids <- function(resids, B,boot_mtd, seed_addition){
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
    }
    T_star[b] = sum(resids * W)
  }
  return(T_star)
}

boot_naive <- function( df_dgp,
                        B,
                        seed_addition,
                        dgp_name="otsu",
                        mu_model="linear") {
  # B <- 100; seed_addition <- 1
  # dgp_name <- "toy"; mu_model = "linear"
  # df_dgp <- gen_one_toy()
  bootstrap_estimates <- numeric(B)

  # X_names<- c("X1","X2")
  # if (mu_model == "linear"){
  #   SL_lib <- "SL.lm"
  # }else if (mu_model == "non-linear"){
  #   SL_lib <-  "SL.randomForest"
  # }
  # # Get the debiased model of the original data
  # SL_fit <- get_SL_fit(df_to_fit=df_dgp,
  #                      X_names = X_names,
  #                      Y_name = "Y",
  #                      SL.library = SL_lib)


  # Perform bootstrap
  for (b in 1:B) {
    # b <- 2
    print(paste0("Naive boot: ", b))
    set.seed(123 + seed_addition + b*13)

    sample_separately=T
    if (sample_separately){
      # Split the data into treated and control groups
      treated <- df_dgp[df_dgp$Z == 1, ]
      control <- df_dgp[df_dgp$Z == 0, ]

      # Resample each group separately
      resampled_treated <- treated[sample(nrow(treated), replace = TRUE), ]
      resampled_control <- control[sample(nrow(control), replace = TRUE), ]

      # Combine the resampled groups
      resampled_data <- rbind(resampled_treated, resampled_control)
    }else{
      resampled_data <-
        df_dgp[sample(nrow(df_dgp), replace = TRUE), ]
    }

    resampled_data$id <- df_dgp$id

    if (dgp_name == "toy"){
      scaling_resampled_data <- resampled_data %>%
        summarize(across(starts_with("X"),
                         function(x) {
                           if (is.numeric(x)) 6 / (max(x) - min(x))
                           else 1000
                         }))
    }else if (dgp_name == "kang"){
      scaling_resampled_data <- resampled_data %>%
        summarize(across(starts_with("X"),
                         function(x) {
                           if (is.numeric(x)) 5 / (max(x) - min(x))
                           else 1000
                         }))
    }

    # Fix the debiasing model
    # estimator <- function(df){
    #   preds_csm <- get_matches(dgp_name=dgp_name,
    #                            df_dgp=df,
    #                            scaling=scaling)
    #   preds_csm[,paste0("hat_mu_0")] <-
    #     SL_pred  <- get_SL_pred(SL_fit=SL_fit,
    #                             df_test=preds_csm,
    #                             X_names=X_names)
    #   ## Construct residuals
    #   # First, construct each \tilde \tau_i
    #   tmp0 <- preds_csm %>%
    #     mutate(Y_bias_corrected = Y - hat_mu_0) %>%
    #     group_by(subclass, Z) %>%
    #     summarize(mn = sum(Y_bias_corrected*weights))
    #   tmp <- tmp0 %>%
    #     group_by(subclass) %>%
    #     summarise(tilde_tau = last(mn)-first(mn))
    #   tilde_tau = tmp$tilde_tau
    #   return(tilde_tau)
    # }

    # Estimate debiasing model independently
    estimator <- function(df){
      obj <- get_matches_and_debiased_residuals(
        dgp_name=dgp_name,
        df_dgp=df,
        scaling=scaling_resampled_data,
        mu_model=mu_model,
        n_split=1)
      return(obj$mean_tilde_tau)
    }
    bootstrap_estimates[b] <-
      estimator(resampled_data)
    print(bootstrap_estimates[b])
  }

  return(bootstrap_estimates)
}





boot_cluster <-function(df_dgp,
                        B,
                        seed_addition,
                        dgp_name="otsu",
                        mu_model="linear"){

  # function(df_dgp,
  #                                 B,
  #                                 seed_addition){
  # B <- 100; seed_addition <- 1
  # dgp_name <- "toy"; mu_model = "linear"
  # df_dgp <- gen_one_toy()
  scaling <- df_dgp %>%
    summarize(across(starts_with("X"),
                     function(x) {
                       if (is.numeric(x)) 6 / (max(x) - min(x))
                       else 1000
                     }))
  obj <- get_matches_and_debiased_residuals(
    dgp_name=dgp_name,
    df_dgp=df_dgp,
    scaling=scaling,
    mu_model="linear",
    n_split=1)
  md <- obj$preds_csm
  pair_ids <- unique(md$subclass)

  #Unit IDs, split by pair membership
  split_inds <-
    split(seq_len(nrow(md)), md$subclass)
  # debiased_units <- md %>%
  #   mutate(Y_bias_corrected = Y - hat_mu_0) %>%
  #   group_by(subclass, Z) %>%
  #   summarize(mn = sum(Y_bias_corrected*weights))%>%
  #   group_by(subclass) %>%
  #   summarise(tilde_tau = last(mn)-first(mn))
  # tilde_tau <- debiased_units$tilde_tau
  # hist(tilde_tau)
  # # Bootstrap process
  # for(i in 1:200) {
  #   sample_data <- sample(tilde_tau,
  #                         replace = TRUE,
  #                         size = length(tilde_tau))
  #   bootstrap_estimates[i] <- mean(sample_data)
  # }
  # hist(bootstrap_estimates)
  # sd(bootstrap_estimates)
  #
  # # Compute bootstrap standard error
  # bootstrap_se <- sd(bootstrap_estimates)



  bootstrap_estimates <- numeric(B)
  # Perform bootstrap
  for (b in 1:B) {
    # b<-1
    # print(b)
    set.seed(123 + seed_addition + b*13)

    n_ids <- length(pair_ids)

    i <- sample(1:n_ids,replace=T)

    ids <- unlist(split_inds[pair_ids[i]])

    resampled_data <- md[ids,]

    estimator <- function(df){
      # Calculate mean_tilde_tau
      tmp0 <- df %>%
        mutate(Y_bias_corrected = Y - hat_mu_0) %>%
        group_by(subclass, Z) %>%
        summarize(mn = sum(Y_bias_corrected*weights))
      tmp <- tmp0 %>%
        group_by(subclass) %>%
        summarise(tilde_tau = last(mn)-first(mn), .groups="drop")
      tilde_tau = tmp$tilde_tau
      mean_tilde_tau <- mean(tilde_tau)
      return(mean_tilde_tau)
    }
    bootstrap_estimates[b] <-
      estimator(resampled_data)
  }

  return(bootstrap_estimates)
}

boot_cluster_for_otsu <- function(df_dgp,
                                  B,
                                  seed_addition){
  # B <- 100; seed_addition <- 1
  # df_dgp <- generate_one_otsu()
  obj <- get_matches_and_debiased_residuals(
    dgp_name="otsu",
    df_dgp=df_dgp,
    scaling=0,
    mu_model="linear",
    n_split=1)
  md <- obj$preds_csm
  pair_ids <- unique(md$subclass)

  #Unit IDs, split by pair membership
  split_inds <-
    split(seq_len(nrow(md)), md$subclass)


  bootstrap_estimates <- numeric(B)
  # Perform bootstrap
  for (b in 1:B) {
    # b<-1
    # print(b)
    set.seed(123 + seed_addition + b*13)

    n_ids <- length(pair_ids)

    i <- sample(1:n_ids,replace=T)

    ids <- unlist(split_inds[pair_ids[i]])

    resampled_data <- md[ids,]

    estimator <- function(df){
      # Calculate mean_tilde_tau
      tmp0 <- df %>%
        mutate(Y_bias_corrected = Y - hat_mu_0) %>%
        group_by(subclass, Z) %>%
        summarize(mn = sum(Y_bias_corrected*weights))
      tmp <- tmp0 %>%
        group_by(subclass) %>%
        summarise(tilde_tau = last(mn)-first(mn), .groups="drop")
      tilde_tau = tmp$tilde_tau
      mean_tilde_tau <- mean(tilde_tau)
      return(mean_tilde_tau)
    }
    bootstrap_estimates[b] <-
      estimator(resampled_data)
  }

  return(bootstrap_estimates)
}


assign_group_label <- function(df_to_split, n_split) {
  df_to_split$group_label <-
    sample( df_to_split$id,
           nrow(df_to_split),
           replace = TRUE)
  return(df_to_split)
}



#### Primary simulation function ####

#' Run a simulation to check inference
#'
#' @param dgp_name Name of the DGP
#' @param att0 True ATT
#' @param I Number of MC runs
#' @param B Number of bootstrap samples per MC run
#' @param mu_model Model for the potential outcomes
#' @param boot_mtd Method for the bootstrap
#' @param n_split Number of splits for sample splitting
#' @param kang_true Whether to use the true ATT from Kang et al. (2021)
#'
#' @return
#' @export
#'
#' @examples
#' boot_otsu_to_test <-
#' boot_CSM(dgp_name="otsu", att0=T,I=100, B=100,mu_model="linear", boot_mtd="Bayesian",n_split=2)
boot_CSM <- function(dgp_name,
                     att0,
                     I=100,
                     B=250,
                     mu_model="linear",
                     boot_mtd="Bayesian",
                     n_split=1,
                     kang_true=FALSE,
                     toy_ctr_dist=0.5,
                     seed = NULL ){

  covered <- CI_lower <- CI_upper <-
    att_true <- att_est <- att_debiased <-
    sd_boot <-
    time_on_matching <-
    time_on_bootstrap <-
    noise <-
    bias <-
    N_T <- N_C <- numeric(I)
  T_star <- numeric(B)
  if ( !is.null(seed) ) {
    set.seed(123)
  }

  for (i in 1:I){
    print(i)
    dgp_obj <-
      get_df_scaling_from_dgp_name(dgp_name=dgp_name,
                                   kang_true=kang_true,
                                   toy_ctr_dist=toy_ctr_dist)
    list2env(dgp_obj,envir = environment())

    time_before_matching <- proc.time()
    matches_and_debiased_residuals<-
      get_matches_and_debiased_residuals(
        dgp_name, df_dgp,
        scaling, mu_model,n_split)
    list2env(matches_and_debiased_residuals,
             envir = environment())
    time_after_matching <- proc.time()
    time_on_matching[i] <-
      (time_after_matching - time_before_matching)[3]

    att_debiased[i] <- mean_tilde_tau

    # Perform bootstrap
    seed_addition = i * 11
    if (boot_mtd == "Bayesian" || boot_mtd == "wild"){
      T_star <-
        boot_by_resids(resids=tilde_tau_resids,
                       B=B,
                       boot_mtd=boot_mtd,
                       seed_addition=seed_addition)

      CI_lower[i] = mean_tilde_tau - quantile(T_star, 0.975)
      CI_upper[i] = mean_tilde_tau - quantile(T_star, 0.025)
      sd_boot[i] = sd(T_star)
    } else if (boot_mtd == "naive"){
      T_star <- boot_naive( df_dgp = df_dgp,
                            B = B,
                            seed_addition = seed_addition,
                            dgp_name=dgp_name,
                            mu_model=mu_model)

      sd_boot[i] = sd(T_star)
      CI_lower[i] = mean_tilde_tau -1.96 *  sd_boot[i]
      CI_upper[i] = mean_tilde_tau + 1.96 * sd_boot[i]
    }else if (boot_mtd == "cluster"){
      T_star <- boot_cluster(df_dgp = df_dgp,
                             B = B,
                             seed_addition = seed_addition,
                             dgp_name=dgp_name,
                             mu_model=mu_model)

      sd_boot[i] = sd(T_star)
      CI_lower[i] = mean_tilde_tau -1.96 *  sd_boot[i]
      CI_upper[i] = mean_tilde_tau + 1.96 * sd_boot[i]
    }else if(boot_mtd == "regression"){
      library(estimatr)
      lm_test_weighted<-
        lm_robust(Y ~ Z,
                  data=preds_csm,
                  weights=weights)
      tmp<-summary(lm_test_weighted)$coefficients
      sd_boot[i] = tmp[2,2]
      CI_lower[i] = tmp[2,5]
      CI_upper[i] = tmp[2,6]
    }else if(boot_mtd == "regression_debiased"){
      library(estimatr)
      preds_csm <- preds_csm %>%
        mutate(Y_bias_corrected = Y - hat_mu_0)
      lm_test_weighted<-
        lm_robust(Y_bias_corrected ~ Z,
                  data=preds_csm,
                  weights=weights)
      tmp<-summary(lm_test_weighted)$coefficients
      sd_boot[i] = tmp[2,2]
      CI_lower[i] = tmp[2,5]
      CI_upper[i] = tmp[2,6]
    }else if (boot_mtd == "A-E"){
      sd_boot[i] = get_se_AE(preds_csm)
      CI_lower[i] = mean_tilde_tau -1.96 *  sd_boot[i]
      CI_upper[i] = mean_tilde_tau + 1.96 * sd_boot[i]
    }

    time_after_bootstrap <- proc.time()
    time_on_bootstrap[i] <-
      (time_after_bootstrap - time_after_matching)[3]


    # Get true ATT for coverage
    if (att0){
      att_true[i] <- att <- 0
    }else{
      att_true[i] <- att <- df_dgp %>%
        filter(Z & !(id %in% attr(preds_csm, "unmatched_units"))) %>%
        summarize(att = mean(Y1-Y0)) %>%
        pull(att)
    }



    covered[i] = (CI_lower[i] < att) & (att < CI_upper[i])




    # Get the CSM estimate
    N_C[i] <- preds_csm %>%
      filter(!(id %in% attr(preds_csm, "unmatched_units"))) %>%
      group_by(Z) %>%
      summarize(sum_weights = sum(weights)) %>%
      slice(1) %>%
      pull(sum_weights)

    N_T[i] <- preds_csm %>%
      filter(!(id %in% attr(preds_csm, "unmatched_units"))) %>%
      group_by(Z) %>%
      summarize(sum_weights = sum(weights)) %>%
      slice(2) %>%
      pull(sum_weights)

    att_est[i] <- preds_csm %>%
      filter(!(id %in% attr(preds_csm, "unmatched_units"))) %>%
      group_by(Z) %>%
      summarize(mn = sum(Y*weights) / sum(weights)) %>%
      summarize(est = last(mn) - first(mn)) %>%
      pull(est)

    noise[i] <- preds_csm %>%
      filter(!(id %in% attr(preds_csm, "unmatched_units"))) %>%
      group_by(Z) %>%
      summarize(mn = sum(noise*weights) / sum(weights)) %>%
      summarize(est = last(mn) - first(mn)) %>%
      pull(est)

    bias[i] <- preds_csm %>%
      filter(!(id %in% attr(preds_csm, "unmatched_units"))) %>%
      group_by(Z) %>%
      summarize(mn = sum(Y0*weights) / sum(weights)) %>%
      summarize(est = last(mn) - first(mn)) %>%
      pull(est)

    print(paste0("att_est is ", signif(att_est[i],3),
                 "; LB is ", signif(CI_lower[i],3),
                 "; UB is ", signif(CI_upper[i],3)))
    print(paste0("Bootstrap s.e. is ", sd_boot[i]))
    print(paste0("Covered is ", covered[i]))

  }
  mean(covered)

  res_save_bayesian_boot <- tibble(id=1:I,
                                   name=dgp_name,
                                   boot_mtd = boot_mtd,
                                   att_true = att_true,
                                   att_est= att_est,
                                   att_est_debiased = att_debiased,
                                   lower=CI_lower,upper=CI_upper,
                                   covered=covered,
                                   sd_boot=sd_boot,
                                   n_split=n_split,
                                   mu_model=mu_model,
                                   time_on_matching=time_on_matching,
                                   time_on_bootstrap=time_on_bootstrap,
                                   noise=noise,
                                   bias=bias,
                                   N_T =N_T,
                                   N_C = N_C)
  return(res_save_bayesian_boot)
}



#' Run a simulation to check inference on the
#' A-E inference method
#'
#' @param dgp_name Name of the DGP
#' @param att0 True ATT
#' @param R Number of MC runs
#' @param B Number of bootstrap samples per MC run
#' @param toy_ctr_dist distance between centers in the to
#'
#' @return
#' @export
#'
#' @examples
sim_inference_CSM_A_E <- function(dgp_name,
                     att0,
                     R=100,
                     B=250,
                     toy_ctr_dist=0.5,
                     seed = NULL){
  # R <- 10; toy_ctr_dist=0.5; dgp_name <- "toy"; att0<-F
  covered <- CI_lower <- CI_upper <-
    att_true <- att_est <- att_debiased <-
    se_AE <-
    time_on_matching <-
    time_on_bootstrap <-
    error <-
    bias <-
    N_T <- N_C <- numeric(R)
  if ( is.null(seed) ) {
    set.seed(123)
  }
  ### Generate one dataset
  df_dgp <-
    gen_one_toy(ctr_dist=toy_ctr_dist) %>%
    mutate(Y0_denoised = Y0 - noise,
           Y1_denoised = Y1 - noise,
           Y_denoised = Y - noise)
  scaling <- 8
  for (i in 1:R){
    # i <- 1
    print(i)

    ## Replace
    df_dgp_i <- df_dgp %>%
      mutate(noise = rnorm(n(), mean=0, sd=0.5)) %>%
      mutate(Y0 = Y0_denoised + noise,
             Y1 = Y1_denoised + noise,
             Y = Y_denoised + noise)
    ### Perform matching
    mtch <- get_cal_matches(
      df_dgp_i,
      metric = "maximum",
      scaling = scaling,
      caliper = 1,
      rad_method = "adaptive",
      est_method = "scm"
    )



    ### Perform inference using the A-E method
    ATT_estimate <- get_ATT_estimate( mtch )
    att_est[i] <- ATT_estimate$ATE
    se_AE[i] = ATT_estimate$SE
    CI_lower[i] = att_est[i] -1.96 *  se_AE[i]
    CI_upper[i] = att_est[i] + 1.96 * se_AE[i]

    treatment_table <- mtch$treatment_table


    ### Obtain necessary things to output
    # Get true ATT
    if (att0){
      att_true[i] <- att <- 0
    }else{
      att_true[i] <- att <-
        df_dgp_i %>%
        filter(Z & (id %in% treatment_table$id)) %>%
        summarize(att = mean(Y1-Y0)) %>%
        pull(att)
    }

    covered[i] = (CI_lower[i] < att) & (att < CI_upper[i])

    # Get error and bias
    full_units <- full_unit_table(mtch, nonzero_weight_only = TRUE )
    error[i] <- full_units %>%
      group_by(Z) %>%
      summarize(mn = sum(noise*weights) / sum(weights)) %>%
      summarize(est = last(mn) - first(mn)) %>%
      pull(est)

    bias[i] <- full_units %>%
      group_by(Z) %>%
      summarize(mn = sum(Y0_denoised*weights) / sum(weights)) %>%
      summarize(est = last(mn) - first(mn)) %>%
      pull(est)

    print(paste0("att_est is ", signif(att_est[i],3),
                 "; LB is ", signif(CI_lower[i],3),
                 "; UB is ", signif(CI_upper[i],3)))
    print(paste0("A-E s.e. is ", se_AE[i]))
    print(paste0("Covered is ", covered[i]))

  }

  res_save_bayesian_boot <-
    tibble(id=1:R,
           name=dgp_name,
           att_true = att_true,
           att_est= att_est,
           lower=CI_lower,
           upper=CI_upper,
           covered=covered,
           se_AE=se_AE,
           error=error,
           bias=bias)
  return(res_save_bayesian_boot)
}



