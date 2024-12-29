
# Code to run the section5.4 simulation



#### Helper functions ####
### Generate Otsu and Rai DGP to debug the bootstrap
m <- function(z) {
  0.4 + 0.25 * sin(8 * z - 5) + 0.4 * exp(-16 * (4 * z - 2.5)^2)
}

gen_df_otsu <- function(N = 1000,K = 2,
                        N1 = NULL,
                        N0 = NULL){
  # output: "id" "X1" "X2" "Z"  "Y0" "Y1" "Y"
  # Right now we choose m to be the 6th specification
  if (!is.null(N1) & !is.null(N0)){
    N <- N1 + N0
  }

  # Generate data
  xi <- runif(N, 0, 1)
  zeta <- matrix(rnorm(N * K), ncol = K)
  X <- xi * apply(abs(zeta), 1, function(x) x / sqrt(sum(x^2)))
  norm_X <- sqrt(colSums(X^2))



  if (!is.null(N1) & !is.null(N0)){
    Z <- c(rep(1, N1), rep(0, N0))
  }else{
    gamma1 <- 0.15
    gamma2 <- 0.7
    P_X <- gamma1 + gamma2 * norm_X

    vi <- runif(N)
    Z <- as.numeric(P_X >= vi)
  }


  # Potential outcomes
  epsilon <- rnorm(N, 0, 0.2)
  Y0 <- m(norm_X) + epsilon
  tau <- 0  # tau is set to 0
  Y1 <- Y0 + tau
  Y <- (1 - Z) * Y0 + Z * Y1

  # Create dataframe
  df <- data.frame(id = 1:N,
                   X1 = X[1, ],
                   X2 = X[2, ],
                   norm_X = norm_X,
                   m_X = m(norm_X),
                   Z = Z,
                   Y0 = Y0,
                   Y1 = Y1,
                   Y = Y)
  return(df)
}

# generate_one_otsu <- function(){
#   # gen_df_otsu(N = 1000,K = 2)
#   # Edit 26 Dec 2023: Realized sample size
#   #   might affect the success of the s.e.
#   # The paper has sample size 100. Thus change
#   #     to it
#   gen_df_otsu(N = 100,K = 2)
# }

get_df_scaling_from_dgp_name <- function(dgp_name,
                                         kang_true = F,
                                         toy_ctr_dist=0.5,
                                         N1 = NULL,
                                         N0 = NULL){
  if (dgp_name == "toy"){
    df_dgp <- gen_one_toy(ctr_dist=toy_ctr_dist)
    scaling <- 8
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
    if (is.null(N1) & is.null(N0)){
      df_dgp <- gen_df_otsu(N = 100,K = 2)
    }else{
      df_dgp <- gen_df_otsu(N1 = N1,
                            N0 = N0,
                            K = 2)
    }

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



debias_residuals <- function(preds_csm, X_names, Y_name, SL_lib, n_split) {
  # Filter controls
  df_controls <- preds_csm %>% filter(Z == 0)

  # Nest data for sample splitting
  df_to_fit_nested <- df_controls %>%
    group_by(group_label) %>%
    nest() %>%
    arrange(group_label)


  fit_and_predict_split <- function(df_to_fit, preds_csm, X_names, Y_name, SL_lib, split_index) {
    # Fit the model
    SL_fit <- get_SL_fit(
      df_to_fit = df_to_fit,
      X_names = X_names,
      Y_name = Y_name,
      SL.library = SL_lib
    )

    # Predict using the fitted model
    preds_csm[, paste0("hat_mu_0_pred_", split_index)] <- get_SL_pred(
      SL_fit = SL_fit,
      df_test = preds_csm,
      X_names = X_names
    )

    return(SL_fit)
  }
  models_n_split <- lapply(seq_len(n_split), function(i) {
    fit_and_predict_split(
      df_to_fit = df_to_fit_nested$data[[i]],
      preds_csm = preds_csm,
      X_names = X_names,
      Y_name = Y_name,
      SL_lib = SL_lib,
      split_index = i
    )
  })

  # Select predictions based on `pred_label`
  preds_csm <- preds_csm %>%
    rowwise() %>%
    mutate(hat_mu_0 = get(paste0("hat_mu_0_pred_", pred_label)))

  # Construct residuals
  tmp0 <- preds_csm %>%
    mutate(Y_bias_corrected = Y - hat_mu_0) %>%
    group_by(subclass, Z) %>%
    summarize(mn = sum(Y_bias_corrected * weights), .groups = "drop")

  tmp <- tmp0 %>%
    group_by(subclass) %>%
    summarise(tilde_tau = last(mn) - first(mn), .groups = "drop")

  tilde_tau <- tmp$tilde_tau
  mean_tilde_tau <- mean(tilde_tau)
  tilde_tau_resids <- tilde_tau - mean_tilde_tau

  return(list(
    preds_csm = preds_csm,
    tilde_tau = tilde_tau,
    mean_tilde_tau = mean_tilde_tau,
    tilde_tau_resids = tilde_tau_resids
  ))
}


get_matches_and_debiased_residuals <- function(
    dgp_name,
    df_dgp,
    scaling,
    mu_model,
    n_split = 1,
    debias = TRUE
) {
  # Step 1: Determine matching type and get matches
  matching_type <- switch(
    dgp_name,
    "toy" = "maximum_fixed_scm",
    "kang" = "maximum_fixed_scm",
    "otsu" = "euclidean_knn",
    stop("Invalid dgp_name. Must be 'toy', 'kang', or 'otsu'.")
  )

  preds_csm <- get_matches(
    matching_type = matching_type,
    df_dgp = df_dgp,
    scaling = scaling
  )

  # Step 2: Assign group labels for cross-fitting
  df_dgp_splitted <- assign_group_label(df_dgp, n_split = n_split) %>%
    mutate(id = as.character(id))

  id_to_group_map <- df_dgp_splitted %>%
    select(id, group_label)

  preds_csm <- left_join(preds_csm, id_to_group_map, by = "id")

  pred_label_map <- get_pred_label_map(n = n_split)
  preds_csm <- left_join(preds_csm, pred_label_map, by = "group_label")

  # Step 3: Debiasing step (optional)
  if (debias) {
    # Determine covariate names and model library
    if (dgp_name == "toy" || dgp_name == "otsu") {
      X_names <- c("X1", "X2")
      SL_lib <- if (mu_model == "linear") "SL.lm" else "SL.randomForest"
    } else if (dgp_name == "kang") {
      if (mu_model == "kang_correct") {
        X_names <- c("V1", "V2", "V3", "V4")
        SL_lib <- "SL.lm"
      } else if (mu_model == "linear") {
        X_names <- c("X1", "X2", "X3", "X4")
        SL_lib <- "SL.lm"
      }
    } else {
      stop("dgp_name must be 'toy' or 'kang'")
    }

    # Perform debiasing
    debiasing_result <- debias_residuals(
      preds_csm = preds_csm,
      X_names = X_names,
      Y_name = "Y",
      SL_lib = SL_lib,
      n_split = n_split
    )

    return(debiasing_result)
  } else {
    return(list(preds_csm = preds_csm))
  }
}

# get_matches_old <- function(dgp_name, df_dgp, scaling){
#   # TODO: WHY should the type of dgp affect how to get matches?  E.g.,
#   # why is otsu different?  It should be a flag for kind of matching,
#   # instead.
#
#   if (dgp_name == "toy" || dgp_name == "kang") {
#     df_dgp_with_matches <- get_cal_matches(
#       df = df_dgp,
#       metric = "maximum",
#       scaling = scaling,
#       rad_method = "fixed",
#       est_method = "scm",
#       return = "all"
#     )
#   } else if (dgp_name == "otsu") {
#     df_dgp_with_matches <- get_cal_matches(
#       df_dgp,
#       covs = starts_with("X"),
#       treatment = "Z",
#       scaling = 1,
#       metric = "euclidean",
#       rad_method = "knn",
#       k = 8
#     )
#     df_dgp_with_matches <- full_unit_table(df_dgp_with_matches,
#                                            nonzero_weight_only = F)
#     # df_dgp_with_matches <- get_NN_matches(df_dgp)
#   } else {
#     stop("dgp_name must be toy, kang, or otsu")
#   }
#   return(df_dgp_with_matches)
# }



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
    sample( 1:n_split,
      # df_to_split$id,
            nrow(df_to_split),
            replace = TRUE)
  return(df_to_split)
}


#' Get Matched Matrix
#'
#' Converts a full matched table into a matrix where each row corresponds to a treated unit
#' and each column contains the IDs of matched control units. Missing values (if any) are filled with `NA`.
#'
#' @param full_matched_table A data frame containing the matched data. Must include the columns:
#' \itemize{
#'   \item{\code{ID}: Unique identifiers for each unit.}
#'   \item{\code{Z}: Indicator for treated (\code{1}) or control (\code{0}) units.}
#'   \item{\code{subclass}: Subclass assignments for matching.}
#' }
#'
#' @return A matrix where:
#' \item{Rows}{Represent treated units (indexed by their IDs).}
#' \item{Columns}{Contain IDs of matched control units for each treated unit. Missing matches are filled with \code{NA}.}
#'
#' @examples
#' # Example full matched table
#' full_matched_table <- data.frame(
#'   ID = 1:6,
#'   Z = c(1, 0, 1, 0, 0, 1),
#'   subclass = c(1, 1, 2, 2, 2, 3)
#' )
#' get_matched_matrix(full_matched_table)
#'
#' @export
get_matched_matrix <- function(full_matched_table) {
  subclasses <- unique(full_matched_table$subclass)

  matched_control_ids <- list()

  for (subclass in subclasses) {
    data_in_subclass <- full_matched_table[full_matched_table$subclass == subclass, ]

    treated_units <- data_in_subclass[data_in_subclass$Z == 1, ]
    control_units <- data_in_subclass[data_in_subclass$Z == 0, ]
    # print(head(treated_units) )
    for (treated_id in treated_units$id) {
      matched_control_ids[[as.character(treated_id)]] <- control_units$id
    }
  }

  max_controls <- max(sapply(matched_control_ids, length))
  matched_matrix <- do.call(rbind, lapply(matched_control_ids, function(ids) {
    c(ids, rep(NA, max_controls - length(ids)))
  }))

  rownames(matched_matrix) <- names(matched_control_ids)

  return(matched_matrix)
}

#' Compute Shared Neighbors
#'
#' Computes the number of shared control neighbors between each pair of treated units
#' based on a given matched matrix.
#'
#' @param matched_matrix A matrix where:
#' \item{Rows}{Represent treated units.}
#' \item{Columns}{Contain IDs of matched control units. Missing values should be \code{NA}.}
#'
#' @return A symmetric matrix of dimensions \code{N1 x N1}, where \code{N1} is the number of treated units.
#' Each entry \code{[i, j]} represents the count of shared control neighbors between treated units \code{i} and \code{j}.
#'
#' @examples
#' # Example matched matrix
#' matched_matrix <- matrix(
#'   c(2, 3, NA, 3, 4, 5, NA, NA, 5, 6, 7, 8),
#'   nrow = 4, byrow = TRUE
#' )
#' compute_shared_neighbors(matched_matrix)
#'
#' @export
compute_shared_neighbors <- function(matched_matrix) {
  N1 <- nrow(matched_matrix)

  shared_neighbors <- matrix(FALSE, nrow = N1, ncol = N1)

  for (i in 1:(N1 - 1)) {
    for (j in (i + 1):N1) {
      shared <- intersect(matched_matrix[i, ], matched_matrix[j, ])
      shared_neighbors[i, j] <- shared_neighbors[j, i] <- length(shared)
    }
  }

  return(shared_neighbors)
}


#' Compute Overlap Statistics
#'
#' Computes summary statistics for the overlap of shared neighbors in a matrix.
#' The function calculates the median, 75th percentile, and maximum for the number of shared treated
#' and control neighbors across rows.
#'
#' @param shared_neighbors A numeric matrix of dimensions \eqn{N1 \times M}, where each entry represents
#' the number of shared neighbors between treated and control units.
#'
#' @return A list containing the following components:
#' \item{avg_shared_controls}{Median number of shared control neighbors across rows.}
#' \item{p75_shared_controls}{75th percentile of the number of shared control neighbors across rows.}
#' \item{max_shared_controls}{Maximum number of shared control neighbors across rows.}
#' \item{avg_shared_treated}{Median number of shared treated neighbors across rows.}
#' \item{p75_shared_treated}{75th percentile of the number of shared treated neighbors across rows.}
#' \item{max_shared_treated}{Maximum number of shared treated neighbors across rows.}
#'
#' @examples
#' # Example matrix with shared neighbors
#' shared_neighbors <- matrix(c(0, 1, 2, 3, 0, 1, 4, 5, 2, 1), nrow = 5, byrow = TRUE)
#' compute_overlap_statistics(shared_neighbors)
#'
#' @export
compute_overlap_statistics <- function(shared_neighbors) {
  N1 <- nrow(shared_neighbors)

  shared_neighbors_binary <- shared_neighbors > 0
  n_shared_treated_vec <- rowSums(shared_neighbors_binary)
  n_shared_controls_vec <- rowSums(shared_neighbors)

  list(
    avg_shared_controls = median(n_shared_controls_vec),
    p75_shared_controls = quantile(n_shared_controls_vec, probs = 0.75),
    max_shared_controls = max(n_shared_controls_vec),
    avg_shared_treated = median(n_shared_treated_vec),
    p75_shared_treated = quantile(n_shared_treated_vec, probs = 0.75),
    max_shared_treated = max(n_shared_treated_vec)
  )
}




#### Primary simulation function ####

#' Run a simulation to check inference
#   Repeat I times of DGP --> match --> getting bootstrapped intervals
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
                     M = 8,
                     seed = NULL,
                     N1 = NULL,
                     N0 = NULL){

  covered <- CI_lower <- CI_upper <-
    att_true <- att_est <- att_debiased <-
    sd_boot <-
    time_on_matching <-
    time_on_bootstrap <-
    noise <-
    bias <-
    N_T <- N_C <-
    avg_shared_controls <-
    p75_shared_controls <-
    avg_shared_treated <-
    p75_shared_treated <- numeric(I)
  T_star <- numeric(B)
  if ( !is.null(seed) ) {
    set.seed(123)
  }

  for (i in 1:I){
    print(i)
    dgp_obj <-
      get_df_scaling_from_dgp_name(
        dgp_name=dgp_name,
        kang_true=kang_true,
        toy_ctr_dist=toy_ctr_dist,
        N1 = N1,
        N0 = N0)
    df_dgp <- dgp_obj$df_dgp
    scaling <- dgp_obj$scaling

    time_before_matching <- proc.time()
    # input: DGP, scaling (scaling does not matter if we do M-NN matching)
    # output: the matched object mtch
    # We want to do a M-NN matching using M=8.
    # Then we want to do the debiasing step
    # Right now the
    matches_and_debiased_residuals <- get_matches_and_debiased_residuals(
      dgp_name = dgp_name,
      df_dgp = df_dgp,
      scaling = scaling,
      mu_model = mu_model,
      n_split = n_split
    )

    # Explicitly unpack the returned list
    preds_csm <- matches_and_debiased_residuals$preds_csm
    tilde_tau <- matches_and_debiased_residuals$tilde_tau
    mean_tilde_tau <- matches_and_debiased_residuals$mean_tilde_tau
    tilde_tau_resids <- matches_and_debiased_residuals$tilde_tau_resids

    time_after_matching <- proc.time()
    time_on_matching[i] <-
      (time_after_matching - time_before_matching)[3]

    att_debiased[i] <- mean_tilde_tau

    # Calculate overlap metrics
    matched_matrix <- get_matched_matrix(preds_csm)
    shared_neighbors <- compute_shared_neighbors(matched_matrix)
    overlap_statistics <- compute_overlap_statistics(shared_neighbors)
    avg_shared_controls[i] <- overlap_statistics[["avg_shared_controls"]]
    p75_shared_controls[i] <- overlap_statistics[["p75_shared_controls"]]
    avg_shared_treated[i] <- overlap_statistics[["avg_shared_treated"]]
    p75_shared_treated[i] <- overlap_statistics[["p75_shared_treated"]]

    # Perform bootstrap
    # Input: tilde_tau_resids: a vector of length n of
    # Output: sd, confint
    seed_addition = i * 11
    if (boot_mtd == "Bayesian" || boot_mtd == "wild" || boot_mtd == "naive-resid"){
      # T_star <-
      #   boot_by_resids(resids=tilde_tau_resids,
      #                  B=B,
      #                  boot_mtd=boot_mtd,
      #                  seed_addition=seed_addition)
      #
      # CI_lower[i] = mean_tilde_tau - quantile(T_star, 0.975)
      # CI_upper[i] = mean_tilde_tau - quantile(T_star, 0.025)
      # sd_boot[i] = sd(T_star)
      boot_ci <- make_bootstrap_ci(boot_mtd)  # or "wild" or "naive-resid"
      results <- boot_ci(
        resids = tilde_tau_resids,
        mean_est = mean_tilde_tau,
        B = B,
        seed_addition = seed_addition
      )

      CI_lower[i] <- results$ci_lower
      CI_upper[i] <- results$ci_upper
      sd_boot[i] <- results$sd
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

  res_save_bayesian_boot <- tibble(
    id=1:I,
    name=dgp_name,
    boot_mtd = boot_mtd,
    att_true = att_true,
    att_est= att_est,
    att_est_debiased = att_debiased,
    lower=CI_lower,
    upper=CI_upper,
    covered=covered,
    sd_boot=sd_boot,
    n_split=n_split,
    mu_model=mu_model,
    time_on_matching=time_on_matching,
    time_on_bootstrap=time_on_bootstrap,
    noise=noise,
    bias=bias,
    N_T =N_T,
    N_C = N_C,
    avg_shared_controls = avg_shared_controls,
    p75_shared_controls = p75_shared_controls,
    avg_shared_treated = avg_shared_treated,
    p75_shared_treated = p75_shared_treated)

  return(res_save_bayesian_boot)
}


### Result reporting ###
create_bootstrap_comparison_plot <-
  function(
    boot_otsu_wild,
    boot_otsu_A_E,
    output_path = here("scripts/boot/figures/ci_comparison_plot.png")) {
  library(tidyverse)

  # Define true values for reference
  true_coverage <- 0.9473
  true_CI_length <- 0.2381

  # Calculate statistics and summarize results
  wild_results <- boot_otsu_wild %>%
    mutate(CI_length = upper - lower) %>%
    summarise(
      coverage = mean(covered),
      avg_CI_length = mean(CI_length)
    ) %>%
    mutate(method = "Wild Bootstrap")

  A_E_results <- boot_otsu_A_E %>%
    mutate(CI_length = upper - lower) %>%
    summarise(
      coverage = mean(covered),
      avg_CI_length = mean(CI_length)
    ) %>%
    mutate(method = "A-E Bootstrap")

  cat("Wild Bootstrap Results:\n")
  cat("Coverage:", wild_results$coverage, "\n")
  cat("Average CI Length:", wild_results$avg_CI_length, "\n\n")

  cat("A-E Bootstrap Results:\n")
  cat("Coverage:", A_E_results$coverage, "\n")
  cat("Average CI Length:", A_E_results$avg_CI_length, "\n\n")

  # Overlapping statistics
  cat("### Overlapping Statistics ###\n")
  cat("Wild Bootstrap:\n")
  cat("Average Shared Treated Units:", mean(boot_otsu_wild$avg_shared_treated), "\n")
  cat("Average Shared Control Units:", mean(boot_otsu_wild$avg_shared_controls), "\n\n")

  cat("A-E Bootstrap:\n")
  cat("Average Shared Treated Units:", mean(boot_otsu_A_E$avg_shared_treated), "\n")
  cat("Average Shared Control Units:", mean(boot_otsu_A_E$avg_shared_controls), "\n\n")


  # Combine results with true values for plotting
  results <- bind_rows(wild_results, A_E_results) %>%
    mutate(
      coverage_diff = coverage - true_coverage,
      CI_length_diff = avg_CI_length - true_CI_length
    )

  # Include true values as a separate row for reference
  true_results <- tibble(
    method = "True Values",
    coverage = true_coverage,
    avg_CI_length = true_CI_length,
    coverage_diff = 0,
    CI_length_diff = 0
  )

  results <- bind_rows(results, true_results)

  # Create a bar plot for comparison
  results_long <- results %>%
    pivot_longer(
      cols = c(coverage, avg_CI_length),
      names_to = "metric",
      values_to = "value"
    )

  plot <- ggplot(results_long, aes(x = method, y = value, fill = metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      title = "Comparison of Wild Bootstrap, A-E Bootstrap, and True Values",
      x = "Method",
      y = "Value",
      fill = "Metric"
    ) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") +
    geom_hline(aes(yintercept = true_coverage), linetype = "dashed", color = "blue") +
    geom_hline(aes(yintercept = true_CI_length), linetype = "dotted", color = "red") +
    annotate("text", x = 2.5, y = true_coverage + 0.02, label = "True Coverage", color = "blue") +
    annotate("text", x = 2.5, y = true_CI_length - 0.02, label = "True CI Length", color = "red")

  # Save the plot
  ggsave(output_path, plot, width = 8, height = 6)

  # Return the plot object for further use or inspection
  return(plot)
}

