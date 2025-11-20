# scripts/sims/canonical_run_comprehensive_simulation.R
# Modular function to run ALL methods (old + new) on any simulation DGP

#' Run comprehensive simulation with all methods
#'
#' @param sim_type One of "acic", "hainmueller", "kang"
#' @param n_iter Number of iterations
#' @param output_dir Directory to save results
#' @param seed_start Starting seed (default 1000)
run_comprehensive_sim <- function(sim_type,
                                  n_iter = 100,
                                  output_dir = NULL,
                                  seed_start = 1000) {

  require(tidyverse)
  require(mvtnorm)
  require(optweight)
  require(dbarts)
  require(tmle)
  require(AIPW)
  require(grf)
  require(twang)
  require(kbal)
  require(tictoc)
  require(CSM)

  # Set up output directory
  if (is.null(output_dir)) {
    output_dir <- here::here("data", "outputs", "sim_comprehensive", sim_type)
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Load helper functions
  source(here::here("scripts/wrappers.R"))
  source(here::here("R/utils.R"))
  source("R/sim_data.R")

  # Set superlearner libraries
  SL.library1 <- c("SL.mean", "SL.lm", "SL.glm")
  SL.library2 <- c("SL.glm", "SL.glmnet", "SL.randomForest", "SL.xgboost")
  SL.library3Q <- c("SL.glm", "tmle.SL.dbarts2", "SL.glmnet")
  SL.library3g <- c("SL.glm", "tmle.SL.dbarts.k.5", "SL.gam")

  cat("\n========================================\n")
  cat(sprintf("Running %s simulation with ALL methods\n", toupper(sim_type)))
  cat(sprintf("Iterations: %d\n", n_iter))
  cat("========================================\n\n")

  # Store all results
  all_results <- list()

  for (i in 1:n_iter) {

    cat("\n========================================\n")
    cat(sprintf("Iteration %d of %d\n", i, n_iter))
    cat("========================================\n")

    # Set seed
    if (sim_type == "acic") {
      # ACIC uses random seeds
      current_seed <- sample(1:100000, 1)
      set.seed(current_seed)
    } else {
      current_seed <- seed_start + i
      set.seed(current_seed)
    }

    iter_start <- Sys.time()

    # ========================================================================
    # GENERATE DATA - Different DGPs for different simulations
    # ========================================================================
    cat("Generating data...\n")
    data_start <- Sys.time()

    if (sim_type == "acic") {
      # ACIC DGP
      n <- 1000
      p <- 10
      model.trt <- "step"
      model.rsp <- "step"

      df <- gen_df_acic(
        model.trt = model.trt,
        root.trt = 0.35,
        overlap.trt = "full",
        model.rsp = model.rsp,
        alignment = 0.75,
        te.hetero = "high",
        random.seed = current_seed,
        n = n,
        p = p
      )

      covs <- df %>%
        select(starts_with("X")) %>%
        colnames()

      zform1 <- as.formula(paste0("Z ~ ", paste0(covs, collapse = "+")))
      zform2 <- as.formula(paste0("Z ~ (", paste0(covs, collapse = "+"), ")^2"))
      form1 <- as.formula(paste0("Y ~ ", paste0(covs, collapse = "+")))
      form2 <- as.formula(paste0("Y ~ (", paste0(covs, collapse = "+"), ")^2"))
      nbins <- 5

    } else if (sim_type == "hainmueller") {
      # Hainmueller DGP
      nc <- 250
      nt <- 50

      df <- gen_df_hain(
        nc = nc,
        nt = nt,
        sigma_e = "n100",
        outcome = "nl2",
        sigma_y = 1,
        ATE = 0
      )

      covs <- c("X1", "X2", "X3", "X4", "X5", "X6")
      zform1 <- as.formula("Z ~ X1+X2+X3+X4+X5+X6")
      zform2 <- as.formula("Z ~ (X1+X2+X3+X4+X5+X6)^2")
      form1 <- as.formula("Y ~ X1+X2+X3+X4+X5+X6")
      form2 <- as.formula("Y ~ (X1+X2+X3+X4+X5+X6)^2")
      nbins <- 5

    } else if (sim_type == "kang") {
      # Kang & Schafer DGP
      n <- 1000

      df <- gen_df_kang(n = n)

      covs <- df %>%
        select(starts_with("X")) %>%
        colnames()

      zform1 <- as.formula(paste0("Z ~ ", paste0(covs, collapse = "+")))
      zform2 <- as.formula(paste0("Z ~ (", paste0(covs, collapse = "+"), ")^2"))
      form1 <- as.formula(paste0("Y ~ ", paste0(covs, collapse = "+")))
      form2 <- as.formula(paste0("Y ~ (", paste0(covs, collapse = "+"), ")^2"))
      nbins <- 5

    } else {
      stop("sim_type must be one of: 'acic', 'hainmueller', 'kang'")
    }

    # Common processing
    dist_scaling <- df %>%
      summarize(across(
        starts_with("X"),
        function(x) {
          if (is.numeric(x)) nbins / (max(x) - min(x))
          else 1000
        }
      ))

    # Get matches
    preds_csm <- get_cal_matches(
      df = df,
      metric = "maximum",
      scaling = dist_scaling,
      rad_method = "fixed",
      est_method = "average",
      k = 25
    )

    preds_cem <- get_cem_matches(
      df = df,
      num_bins = nbins,
      est_method = "average",
      return = "all"
    )

    # Record infeasible units
    ninf <- length(attr(preds_csm, "unmatched_units"))
    ninf_cem <- sum(df$Z) - length(attr(preds_cem, "feasible_units"))

    # Filter
    df <- df %>%
      filter(!id %in% attr(preds_csm, "unmatched_units"))

    data_time <- as.numeric(difftime(Sys.time(), data_start, units = "secs"))
    cat(sprintf("  Data generation: %.2f seconds\n", data_time))

    # ========================================================================
    # COMPUTE TRUE ATT
    # ========================================================================
    if (sim_type == "acic") {
      true_ATT <- df %>%
        filter(Z) %>%
        summarize(att = mean(Y1 - Y0)) %>%
        pull(att)
    } else {
      true_ATT <- 0
    }

    cat(sprintf("  True ATT: %.3f\n", true_ATT))
    cat(sprintf("  Sample size: n=%d, n_treated=%d, n_control=%d\n",
                nrow(df), sum(df$Z), sum(!df$Z)))
    cat(sprintf("  Infeasible: CSM=%d, CEM=%d\n", ninf, ninf_cem))

    # ========================================================================
    # COMPUTE ALL ESTIMATES
    # ========================================================================

    # Simple difference
    cat("\nSimple difference...\n")
    t_start <- Sys.time()
    att_diff <- get_att_diff(df)
    t_diff <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat(sprintf("  Estimate: %.3f, Time: %.2fs\n", att_diff, t_diff))

    # Balance methods
    cat("\nBalance methods...\n")
    t_start <- Sys.time()
    att_bal1 <- get_att_bal(df, zform1,
                            rep(0.01, length(covs)))
    t_bal1 <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat(sprintf("  bal1: %.3f, Time: %.2fs\n", att_bal1, t_bal1))

    t_start <- Sys.time()
    att_bal2 <- get_att_bal(df, zform2,
                            rep(0.1, length(covs) + choose(length(covs), 2)))
    t_bal2 <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat(sprintf("  bal2: %.3f, Time: %.2fs\n", att_bal2, t_bal2))

    # Outcome regression
    cat("\nOutcome regression...\n")
    t_start <- Sys.time()
    att_or_lm <- get_att_or_lm(df, form = form2)
    t_or_lm <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat(sprintf("  or_lm: %.3f, Time: %.2fs\n", att_or_lm, t_or_lm))

    t_start <- Sys.time()
    att_or_bart <- get_att_or_bart(df, covs = covs)
    t_or_bart <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat(sprintf("  or_bart: %.3f, Time: %.2fs\n", att_or_bart, t_or_bart))

    # Propensity score
    cat("\nPropensity score...\n")
    t_start <- Sys.time()
    att_ps_lm <- get_att_ps_lm(df, zform2)
    t_ps_lm <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat(sprintf("  ps_lm: %.3f, Time: %.2fs\n", att_ps_lm, t_ps_lm))

    t_start <- Sys.time()
    att_ps_bart <- get_att_ps_bart(df, covs = covs)
    t_ps_bart <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat(sprintf("  ps_bart: %.3f, Time: %.2fs\n", att_ps_bart, t_ps_bart))

    # Matching methods
    cat("\nMatching methods...\n")
    t_start <- Sys.time()
    att_csm_scm <- get_att_csm(df, scaling = dist_scaling,
                               est_method = "scm", rad_method = "fixed")
    t_csm_scm <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat(sprintf("  csm_scm: %.3f, Time: %.2fs\n", att_csm_scm, t_csm_scm))

    t_start <- Sys.time()
    att_csm_avg <- get_att_csm(df, scaling = dist_scaling,
                               est_method = "average", rad_method = "fixed")
    t_csm_avg <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat(sprintf("  csm_avg: %.3f, Time: %.2fs\n", att_csm_avg, t_csm_avg))

    t_start <- Sys.time()
    att_cem_scm <- get_att_cem(df, num_bins = nbins,
                               est_method = "scm", estimand = "CEM-ATT")
    t_cem_scm <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat(sprintf("  cem_scm: %.3f, Time: %.2fs\n", att_cem_scm, t_cem_scm))

    t_start <- Sys.time()
    att_cem_avg <- get_att_cem(df, num_bins = nbins,
                               est_method = "average", estimand = "CEM-ATT")
    t_cem_avg <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat(sprintf("  cem_avg: %.3f, Time: %.2fs\n", att_cem_avg, t_cem_avg))

    t_start <- Sys.time()
    att_onenn <- get_att_1nn(df, scaling = dist_scaling)
    t_onenn <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat(sprintf("  onenn: %.3f, Time: %.2fs\n", att_onenn, t_onenn))

    # Doubly robust
    cat("\nDoubly robust methods...\n")
    t_start <- Sys.time()
    att_tmle1 <- get_att_tmle(df, covs = covs,
                              Q.SL.library = SL.library1,
                              g.SL.library = SL.library1)
    t_tmle1 <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat(sprintf("  tmle1: %.3f, Time: %.2fs\n", att_tmle1, t_tmle1))

    t_start <- Sys.time()
    att_aipw1 <- get_att_aipw(df, covs = covs,
                              Q.SL.library = SL.library1,
                              g.SL.library = SL.library1)
    t_aipw1 <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat(sprintf("  aipw1: %.3f, Time: %.2fs\n", att_aipw1, t_aipw1))

    t_start <- Sys.time()
    att_tmle2 <- get_att_tmle(df, covs = covs,
                              Q.SL.library = SL.library3Q,
                              g.SL.library = SL.library3g)
    t_tmle2 <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat(sprintf("  tmle2: %.3f, Time: %.2fs\n", att_tmle2, t_tmle2))

    t_start <- Sys.time()
    att_aipw2 <- get_att_aipw(df, covs = covs,
                              Q.SL.library = SL.library2,
                              g.SL.library = SL.library2)
    t_aipw2 <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat(sprintf("  aipw2: %.3f, Time: %.2fs\n", att_aipw2, t_aipw2))

    # NEW METHODS
    cat("\nNEW METHODS...\n")

    t_start <- Sys.time()
    att_cf <- tryCatch(
      {
        result <- get_att_causal_forest(df, covs = covs)
        cat(sprintf("  causal_forest: %.3f", result))
        result
      },
      error = function(e) {
        cat(sprintf("  causal_forest: ERROR - %s", e$message))
        return(NA_real_)
      }
    )
    t_cf <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat(sprintf(", Time: %.2fs\n", t_cf))

    t_start <- Sys.time()
    att_twang <- tryCatch(
      {
        result <- get_att_twang(df, form = zform1)
        cat(sprintf("  twang: %.3f", result))
        result
      },
      error = function(e) {
        cat(sprintf("  twang: ERROR - %s", e$message))
        return(NA_real_)
      }
    )
    t_twang <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat(sprintf(", Time: %.2fs\n", t_twang))

    t_start <- Sys.time()
    att_kbal <- tryCatch(
      {
        result <- get_att_kbal(df, covs = covs)
        cat(sprintf("  kbal: %.3f", result))
        result
      },
      error = function(e) {
        cat(sprintf("  kbal: ERROR - %s", e$message))
        return(NA_real_)
      }
    )
    t_kbal <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    cat(sprintf(", Time: %.2fs\n", t_kbal))

    # ========================================================================
    # CALCULATE TOTAL TIMING
    # ========================================================================
    total_time <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))

    cat("\n----------------------------------------\n")
    cat(sprintf("Iteration %d complete! Total time: %.2f seconds\n", i, total_time))
    cat("----------------------------------------\n")

    # ========================================================================
    # ASSEMBLE RESULTS
    # ========================================================================
    res <- tibble(
      runid = i,
      seed = current_seed,
      sim_type = sim_type,
      elapsed_time_secs = total_time,
      ninf = ninf,
      ninf_cem = ninf_cem,
      true_ATT = true_ATT,

      # Old methods
      diff = att_diff,
      bal1 = att_bal1,
      bal2 = att_bal2,
      or_lm = att_or_lm,
      or_bart = att_or_bart,
      ps_lm = att_ps_lm,
      ps_bart = att_ps_bart,
      csm_scm = att_csm_scm,
      csm_avg = att_csm_avg,
      cem_scm = att_cem_scm,
      cem_avg = att_cem_avg,
      onenn = att_onenn,
      tmle1 = att_tmle1,
      aipw1 = att_aipw1,
      tmle2 = att_tmle2,
      aipw2 = att_aipw2,

      # New methods
      causal_forest = att_cf,
      twang = att_twang,
      kbal = att_kbal
    )

    # Save iteration results
    FNAME_ITER <- file.path(output_dir, sprintf("run_%03d.csv", i))
    write_csv(res, FNAME_ITER)

    all_results[[i]] <- res
  }

  # ========================================================================
  # COMBINE AND SUMMARIZE
  # ========================================================================
  cat("\n========================================\n")
  cat("Combining results...\n")
  cat("========================================\n")

  combined_results <- bind_rows(all_results) %>%
    write_csv(file.path(dirname(output_dir), paste0(sim_type, "_comprehensive.csv")))

  cat("\n=== Simulation Complete ===\n")
  cat(sprintf("Simulation type: %s\n", sim_type))
  cat(sprintf("Total iterations: %d\n", n_iter))
  cat(sprintf("Total time: %.2f minutes\n", sum(combined_results$elapsed_time_secs)/60))

  cat("\nMethod performance summary (RMSE):\n")
  performance_summary <- combined_results %>%
    pivot_longer(diff:kbal, names_to = "method", values_to = "estimate") %>%
    group_by(method) %>%
    summarize(
      rmse = sqrt(mean((estimate - true_ATT)^2, na.rm = TRUE)),
      bias = mean(estimate - true_ATT, na.rm = TRUE),
      n_na = sum(is.na(estimate))
    ) %>%
    arrange(rmse)

  print(performance_summary, n = Inf)

  return(combined_results)
}
