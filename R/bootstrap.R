# Function: boot_bayesian
# input: 
#   dgp_name: "toy" or "kang",
#   att0 = T or F: T means set true ATT to 0; F means get the true ATT from the 
#   I = 100: number of MC runs
#   B = 250: number of bootstrap samples per MC run
#   mu_model="linear": string 
#   
# output: A tibble
# tibble(id=1:I,
#        name=dgp_name,
#        boot_mtd = "Bayesian",
#        att_true = att_true,
#        att_est= att_est,
#        att_est_debiased = att_debiased,
#        lower=CI_lower,upper=CI_upper, 
#        covered=covered,
#        sd_boot=sd_boot)
boot_bayesian <- function(dgp_name, 
                          att0,
                          I=100,
                          B=250,
                          mu_model="linear"){
  
  covered <- CI_lower <- CI_upper <- 
    att_true <- att_est <- att_debiased<-numeric(I)
  sd_boot <- numeric(I)
  T_star <- T_star_naive <- numeric(B)
  set.seed(123)
  for (i in 1:I){
    # i <- 1
    print(i)
    dgp_obj <- get_df_scaling_from_dgp_name(dgp_name=dgp_name)
    df_dgp <- dgp_obj$df_dgp
    dist_scaling <- dgp_obj$dist_scaling
    
    ## Obtain the debiased residual
    
    
    # Get the matches 
    preds_csm <- get_cal_matches(
      df = df_dgp,
      metric = "maximum",
      dist_scaling = dist_scaling,
      cal_method = "fixed",
      est_method = "scm",
      # return = "agg_co_units",
      return = "all",
      knn = 25)   
    
    
    if (dgp_name == "toy") {
      # Fit the linear model with covariates X1 and X2 for cases where Z==0
      lm_0 <- lm(Y ~ X1 + X2, data = dplyr::filter(df_dgp, Z == 0))
      
      # Predict using the fitted model and update preds_csm for toy
      preds_csm$hat_mu_0 <- predict(lm_0, newdata = data.frame(X1 = preds_csm$X1,
                                                               X2 = preds_csm$X2))
    } else if (dgp_name == "kang") {
      if (mu_model == "linear"){
        
        # Fit the linear model with covariates X1, X2, X3, and X4 for cases where Z==0
        lm_0 <- lm(Y ~ X1 + X2 + X3 + X4, data = dplyr::filter(df_dgp, Z == 0))
        
        # Predict using the fitted model and update preds_csm for kang
        preds_csm$hat_mu_0 <- predict(lm_0, newdata = data.frame(X1 = preds_csm$X1,
                                                                 X2 = preds_csm$X2,
                                                                 X3 = preds_csm$X3,
                                                                 X4 = preds_csm$X4))
      }else if (mu_model == "nonlinear"){
        # Define the library of algorithms
        SL.library2 <- c("SL.glm", 
                         "SL.glmnet",
                         "SL.randomForest", "SL.xgboost")
        
        # Fit the SuperLearner model
        # Assuming df_dgp is the dataset and Z, X1, X2, X3, X4 are columns in df_dgp
        super_learner_model <- SuperLearner(Y = df_dgp$Y[df_dgp$Z == 0],
                                            X = df_dgp[df_dgp$Z == 0, c("X1", "X2", "X3", "X4")],
                                            SL.library = SL.library2,
                                            method = "method.NNLS")
        
        # Predict using the fitted model
        preds_csm$hat_mu_0 <- predict(super_learner_model, 
                                      newdata = preds_csm[, c("X1", "X2", "X3", "X4")],
                                      type = "response")$pred
      }else if (mu_model == "kang_correct"){
        
        # Fit the linear model with covariates X1, X2, X3, and X4 for cases where Z==0
        lm_0 <- lm(Y ~ V1 + V2 + V3 + V4, data = dplyr::filter(df_dgp, Z == 0))
        
        # Predict using the fitted model and update preds_csm for kang
        preds_csm$hat_mu_0 <- predict(lm_0, newdata = data.frame(V1 = preds_csm$V1,
                                                                 V2 = preds_csm$V2,
                                                                 V3 = preds_csm$V3,
                                                                 V4 = preds_csm$V4))
      }
    } 
    else {
      stop("dgp_name must be toy or kang")
    }
    att_est[i] <- get_att_ests(preds_csm)
    
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
    att_debiased[i] <- mean_tilde_tau <- mean(tilde_tau)
    tilde_tau_resids <- tilde_tau - mean_tilde_tau
    
    # Perform Bayesian bootstrap
    for (b in 1:B){
      # W = gtools::rdirichlet(1, alpha=rep(1,100))
      n1 <- length(tilde_tau_resids)
      W = gtools::rdirichlet(1, alpha=rep(1,n1))
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
    
    
    print(paste0("ATT is ", signif(att,3),
                 "; LB is ", signif(CI_lower[i],3),
                 "; UB is ", signif(CI_upper[i],3)))
    covered[i] = (CI_lower[i] < att) & (att < CI_upper[i])
    
    print(paste0("Covered is ", covered[i]))
    
  }
  mean(covered) 
  
  res_save_bayesian_boot <- tibble(id=1:I,
                                   name=dgp_name,
                                   boot_mtd = "Bayesian",
                                   att_true = att_true,
                                   att_est= att_est,
                                   att_est_debiased = att_debiased,
                                   lower=CI_lower,upper=CI_upper, 
                                   covered=covered,
                                   sd_boot=sd_boot)
  return(res_save_bayesian_boot)
}
