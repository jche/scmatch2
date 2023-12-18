# Function: Obtain the debiased residuals 
# input: 
#   df_dgp, dist_scaling, mu_model
# Output:
#   preds_csm: matched df
#   tilde_tau,
#   mean_tilde_tau,
#   tilde_tau_resids

get_SL_fit <- function(df_to_fit,
                       X_names,
                       Y_name,
                       SL.library){
  SuperLearner(Y = df_to_fit[,Y_name,drop=T],
               X = df_to_fit[,X_names],
               SL.library = SL.library)
}

get_SL_pred <- 
  function(SL_fit, df_pred, X_names){
    SL_pred_obj <- predict(SL_fit,
            newdata = df_pred[,X_names])
    return(SL_pred_obj$library.predict)
  }


split_data <- function(df_to_split, n_split) {
  df_to_split$group_label <- 
    sample(1:n_split, 
           nrow(df_to_split), 
           replace = TRUE)
  return(df_to_split)
}

get_matches_and_debiased_residuals <- 
  function(dgp_name,
           df_dgp, 
           dist_scaling, 
           mu_model,
           n_split=1){
    preds_csm <- get_cal_matches(
      df = df_dgp,
      metric = "maximum",
      dist_scaling = dist_scaling,
      cal_method = "fixed",
      est_method = "scm",
      return = "all",
      knn = 25)   
    
    # 2. Get the debiased models. 
    #   Start with n_split == 1. Test on the toy data
    # 2.1 Split preds_csm into n_split pieces. 
    
    # 2.2 Train both model for each split
    #  2.2.1 Assign the right training set, 
    #       Specify the names X variables and 
    #         Y variable
    #   Create nested dataset over 
    
    #  2.2.2 (done) Given the data, X_names, Y_name, 
    #       return a trained model
    #  2.2.3 Attach the trained model correctly 
    # 2.3 Get the fitted value from the other split
    #   Create a mapping: two columns: 
    #     group_label of the data
    #     group_label of the prediction model to be used
    
    
    df_controls <- df_dgp[df_dgp$Z == 0,]
    
    if (dgp_name == "toy") {
      X_names<- c("X1","X2")
      SL_lib <- "SL.lm"
      
    } else if (dgp_name == "kang") {
      if (mu_model == "kang_correct"){
        X_names<- c("V1","V2","V3","V4")
        SL_lib <- "SL.lm"
      }
    } 
    else {
      stop("dgp_name must be toy or kang")
    }
    # Model fitting 
    SL_fit_lm <- get_SL_fit(df_to_fit=df_controls,
                            X_names = X_names,
                            Y_name = "Y",
                            SL.library = SL_lib)
    
    # get predictions:
    preds_csm$hat_mu_0 <-
      SL_pred  <- get_SL_pred(SL_fit=SL_fit_lm,
                           df_pred=preds_csm,
                           X_names=X_names)
    
    ## Construct residuals
    # First, construct each \tilde \tau_i
    tmp0 <- preds_csm %>% 
      mutate(Y_bias_corrected = Y - hat_mu_0) %>%
      group_by(subclass, Z) %>%
      summarize(mn = sum(Y_bias_corrected*weights)) 
    tmp <- tmp0 %>%
      group_by(subclass) %>%
      summarise(tilde_tau = last(mn)-first(mn)) 
    tilde_tau = tmp$tilde_tau
    
    # Second, residual = \tilde \tau_i - mean(\tilde \tau_i)
    mean_tilde_tau <- mean(tilde_tau)
    tilde_tau_resids <- tilde_tau - mean_tilde_tau
    return(list(preds_csm=preds_csm,
        tilde_tau=tilde_tau,
        mean_tilde_tau=mean_tilde_tau,
        tilde_tau_resids=tilde_tau_resids))
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
                     mu_model="linear"){
  boot_CSM(dgp_name, 
           att0,
           I,
           B,
           mu_model,
           boot_mtd="Bayesian")
}

boot_wild <- function(dgp_name, 
                          att0,
                          I=100,
                          B=250,
                          mu_model="linear"){
  boot_CSM(dgp_name, 
           att0,
           I,
           B,
           mu_model,
           boot_mtd="wild")
}

boot_CSM <- function(dgp_name, 
                          att0,
                          I=100,
                          B=250,
                          mu_model="linear",
                          boot_mtd="Bayesian"){
  
  covered <- CI_lower <- CI_upper <- 
    att_true <- att_est <- att_debiased <-
    sd_boot <-numeric(I)
  T_star <- numeric(B)
  set.seed(123)
  for (i in 1:I){
    print(i)
    dgp_obj <- get_df_scaling_from_dgp_name(dgp_name=dgp_name)
    list2env(dgp_obj,envir = environment())
    
    matches_and_debiased_residuals<-
      get_matches_and_debiased_residuals(
         dgp_name, df_dgp, 
        dist_scaling, mu_model)
    list2env(matches_and_debiased_residuals, 
             envir = environment())
    
    att_debiased[i] <- mean_tilde_tau
    
    # Perform Bayesian bootstrap
    for (b in 1:B){
      set.seed(123 + i * 11 + b*13)
      n1 <- length(tilde_tau_resids)
      # The implemented W is W(in the paper) / sqrt(n)
      if (boot_mtd=="Bayesian"){
        W = gtools::rdirichlet(1, alpha=rep(1,n1))
      }else if (boot_mtd=="wild"){
        W = sample(
          c( -(sqrt(5)-1)/2, (sqrt(5)+1)/2 ),
          prob = c( (sqrt(5)+1)/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5)) ),
          replace = T, size = n1) / n1 
      }
      T_star[b] = sum(tilde_tau_resids * W)
      if (length(W) != length(tilde_tau_resids)){
        print(paste0("Length of W is ", length(W),". Length of resid is ",length(tilde_tau_resids)))
      }
      
    }
    
    CI_lower[i] = mean_tilde_tau - quantile(T_star, 0.975)
    CI_upper[i] = mean_tilde_tau - quantile(T_star, 0.025)
    sd_boot[i] = sd(T_star)
    
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
    
    print(paste0("ATT is ", signif(att,3),
                 "; LB is ", signif(CI_lower[i],3),
                 "; UB is ", signif(CI_upper[i],3)))
    print(paste0("Covered is ", covered[i]))
    
    # Get the CSM estimate
    att_est[i] <- get_att_ests(preds_csm)
    
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
                                   sd_boot=sd_boot)
  return(res_save_bayesian_boot)
}
