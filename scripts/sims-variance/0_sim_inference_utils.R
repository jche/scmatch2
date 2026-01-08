# scripts/sims-variance/0_sim_inference_utils.R

devtools::load_all()

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(here)
  library(mvtnorm)
})

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (is.atomic(x) && all(is.na(x)))) y else x
}

load_libs <- function(libs) {
  for (lib in libs) {
    if (!requireNamespace(lib, quietly = TRUE)) {
      warning(paste("Package", lib, "not found. Please install it."), call. = FALSE)
      return(FALSE)
    }
    suppressPackageStartupMessages(library(lib, character.only = TRUE))
  }
  TRUE
}

save_res_to_csv <- function(curr_res, FNAME) {
  if (file.exists(FNAME)) {
    readr::write_csv(curr_res, FNAME, append = TRUE)
  } else {
    readr::write_csv(curr_res, FNAME)
  }
}

get_sim_paths <- function() {
  output_dir_base <- here::here("data/outputs/sims-variance")
  dir.create(output_dir_base, showWarnings = FALSE, recursive = TRUE)
  list(
    output_dir_base = output_dir_base,
    individual_dir = file.path(output_dir_base, "individual"),
    summary_csv = file.path(output_dir_base, "summary_table.csv"),
    combined_csv = file.path(output_dir_base, "combined_results.csv")
  )
}

# Collect rds results into one CSV (optional utility)
collect_results_to_csv <- function(combined_csv = here::here("toy-sim/data/combined_results.csv")) {
  paths <- get_sim_paths()
  individual_dir <- paths$individual_dir
  dir.create(individual_dir, showWarnings = FALSE, recursive = TRUE)

  result_files <- list.files(
    path = individual_dir,
    pattern = "^toy_match_infer_iter_\\d+\\.rds$",
    full.names = TRUE
  )

  if (length(result_files) == 0) {
    stop("No result files found in: ", individual_dir, call. = FALSE)
  }

  combined_results <- tibble()
  failed_files <- character()

  for (fp in result_files) {
    res <- tryCatch(readRDS(fp), error = function(e) NULL)
    if (is.null(res)) {
      failed_files <- c(failed_files, basename(fp))
    } else {
      combined_results <- bind_rows(combined_results, res)
    }
  }

  if (nrow(combined_results) == 0) {
    stop("No valid results were found.", call. = FALSE)
  }

  readr::write_csv(combined_results, combined_csv)
  invisible(list(n_rows = nrow(combined_results), failed = failed_files))
}

# -------------------------------------------------------------------
# Core: Generate toy df (choose fixed vs random mixture for treated/controls)
# -------------------------------------------------------------------

make_csm_toy_df <- function(
    nc = 500, nt = 100,
    f0_sd = 0.5,
    prop_nc_unif = 1/3,
    ctr_dist = 0.5,
    seed = NULL,
    mixture_fixed = FALSE
) {
  if (!is.null(seed)) set.seed(seed)

  # --- original (x,y) versions ---
  f0_fun_xy <- function(x, y) {
    mvtnorm::dmvnorm(
      cbind(x, y),
      mean  = c(0.5, 0.5),
      sigma = matrix(c(1, 0.8, 0.8, 1), 2, 2)
    ) * 20
  }

  tx_effect_fun_xy <- function(x, y) { 3 * x + 3 * y }

  # --- wrappers for gen_df_adv_k (matrix input) ---
  f0_fun_mat <- function(X) {
    X <- as.matrix(X)
    f0_fun_xy(X[, 1], X[, 2])
  }

  tx_effect_fun_mat <- function(X) {
    X <- as.matrix(X)
    tx_effect_fun_xy(X[, 1], X[, 2])
  }

  df <- if (isTRUE(mixture_fixed)) {
    gen_df_adv(
      nc = nc, nt = nt, f0_sd = f0_sd,
      tx_effect_fun = tx_effect_fun_xy,
      f0_fun = f0_fun_xy,
      ctr_dist = ctr_dist,
      prop_nc_unif = prop_nc_unif
    )
  } else {
    gen_df_adv_k(
      nc = nc, nt = nt, k = 2,
      f0_sd = f0_sd,
      tx_effect_fun = tx_effect_fun_mat,
      f0_fun = f0_fun_mat,
      ctr_dist = ctr_dist,
      prop_nc_unif = prop_nc_unif
    )
  }

  df %>% mutate(Z = as.integer(Z))
}


# -------------------------------------------------------------------
# Scaling for toy: follow your snippet
# dist_scaling: (2*nbins)/(range) for numeric covs; nbins=6 for toy
# -------------------------------------------------------------------

compute_toy_scaling <- function(df, nbins = 6) {
  df %>%
    summarise(across(
      starts_with("X"),
      function(x) {
        if (is.numeric(x)) (2 * nbins) / (max(x) - min(x)) else 1000
      }
    ))
}

# -------------------------------------------------------------------
# One iteration: generate data, match, infer
# -------------------------------------------------------------------

toy_match_infer <- function(
    i,
    overlap_label,
    error_label = "homoskedastic",
    prop_nc_unif,
    nc = 500,
    nt = 100,
    f0_sd = 0.5,
    ctr_dist = 0.5,
    nbins = 6,
    include_bootstrap = TRUE,
    boot_mtd = "wild",
    B = 250,
    seed_addition = 11,
    grid_id = 1,
    verbose = FALSE
) {

  # seed for this iteration (already computed in sim_master; pass in seed_addition)
  set.seed(seed_addition)

  # --- generate df (CSM toy spec) ---
  df <- make_csm_toy_df(
    nc = nc, nt = nt,
    f0_sd = f0_sd,
    prop_nc_unif = prop_nc_unif,
    ctr_dist = ctr_dist,
    seed = seed_addition
  )

  covs <- c("X1", "X2")
  scaling <- compute_toy_scaling(df, nbins = nbins)

  # --- matching: mimic your snippet intent ---
  # adaptive caliper w/ k=5 (shrink/expand) and est_method "scm"
  mtch <- tryCatch({
    get_cal_matches(
      data = df,
      Z ~ X1 + X2,
      rad_method = "adaptive",
      scaling = scaling,
      k = 5,
      warn = FALSE,
      est_method = "scm"
    )
  }, error = function(e) {
    if (verbose) message("Matching failed iter ", i, ": ", e$message)
    NULL
  })

  if (is.null(mtch)) {
    message("DEBUG: mtch is NULL (matching failed).")
  } else {
    message("DEBUG: matching succeeded. class(mtch) = ", paste(class(mtch), collapse=", "))
    message("DEBUG: params(mtch)$treatment = ", tryCatch(params(mtch)$treatment, error=function(e) "ERR"))
  }

  rt <- tryCatch(result_table(mtch), error=function(e) NULL)
  rt_sc <- tryCatch(result_table(mtch, "sc_units"), error=function(e) NULL)

  message("DEBUG result_table(mtch) names: ", paste(names(rt), collapse=", "))
  if (!is.null(rt_sc)) message("DEBUG result_table(mtch,'sc_units') names: ", paste(names(rt_sc), collapse=", "))

  # Sanity checks
  if (!is.null(rt)) {
    message("DEBUG nrow(rt)=", nrow(rt))
    message("DEBUG has weights? ", "weights" %in% names(rt))
    message("DEBUG has subclass? ", "subclass" %in% names(rt))
    message("DEBUG has id? ", "id" %in% names(rt))
    message("DEBUG has Y? ", "Y" %in% names(rt))
    message("DEBUG has Z? ", "Z" %in% names(rt))
  }



  # overlap stats (optional; keep light)
  overlap_stats <- list()
  if (!is.null(mtch)) {
    overlap_stats <- tryCatch({
      overlap_obj <- calculate_overlap_statistics_from_match_object(mtch)
      cr <- overlap_obj$control_reuse
      sp <- overlap_obj$shared_per_treated
      list(
        mean_reuse = cr$mean_reuse,
        median_reuse = cr$median_reuse,
        mean_shared_per_treated = sp$mean_shared,
        median_shared_per_treated = sp$median_shared
      )
    }, error = function(e) list())
  }

  # --- inference methods ---
  results_list <- list()

  get_one <- function(method_name, expr) {
    out <- tryCatch(
      expr,
      error = function(e) {
        message("\n==============================")
        message("DEBUG get_ATT_estimate failed for method: ", method_name)
        message("ERROR: ", conditionMessage(e))
        message("TRACEBACK:\n", paste(capture.output(traceback()), collapse = "\n"))
        message("==============================\n")
        return(structure(list(.error = e), class = "att_error"))
      }
    )
    if (inherits(out, "att_error")) return(NULL)
    out <- as_tibble(out)
    out$inference_method <- method_name
    out
  }


  att_pooled <- if (!is.null(mtch)) get_one("pooled", get_ATT_estimate(mtch, variance_method = "pooled")) else NULL
  att_pooled_het <- if (!is.null(mtch)) get_one("pooled_het", get_ATT_estimate(mtch, variance_method = "pooled_het")) else NULL
  att_ai06 <- if (!is.null(mtch)) {
    get_one("ai06", get_ATT_estimate(mtch, variance_method = "ai06", df = df, M = 5, covs = covs))
  } else NULL

  att_boot <- if (include_bootstrap && !is.null(mtch)) {
    get_one("bootstrap", get_ATT_estimate(mtch, variance_method = "bootstrap", boot_mtd = boot_mtd, B = B, seed_addition = seed_addition))
  } else NULL

  extract_tbl <- function(att_tbl, method_name) {
    if (is.null(att_tbl) || nrow(att_tbl) == 0) {
      return(tibble(
        runID = i,
        inference_method = method_name,
        att_est = NA_real_, SE = NA_real_,
        CI_lower = NA_real_, CI_upper = NA_real_,
        N_T = NA_integer_, ESS_C = NA_real_,
        V = NA_real_, V_E = NA_real_, V_P = NA_real_, sigma_hat = NA_real_,
        t = NA_real_
      ))
    }
    tibble(
      runID = i,
      inference_method = method_name,
      att_est = att_tbl$ATT %||% NA_real_,
      SE = att_tbl$SE %||% NA_real_,
      CI_lower = (att_tbl$ATT %||% NA_real_) - 1.96 * (att_tbl$SE %||% NA_real_),
      CI_upper = (att_tbl$ATT %||% NA_real_) + 1.96 * (att_tbl$SE %||% NA_real_),
      N_T = att_tbl$N_T %||% NA_integer_,
      ESS_C = att_tbl$ESS_C %||% NA_real_,
      V = att_tbl$V %||% NA_real_,
      V_E = att_tbl$V_E %||% NA_real_,
      V_P = att_tbl$V_P %||% NA_real_,
      sigma_hat = att_tbl$sigma_hat %||% NA_real_,
      t = att_tbl$t %||% ((att_tbl$ATT %||% NA_real_) / (att_tbl$SE %||% NA_real_))
    )
  }

  results_list[["pooled"]] <- extract_tbl(att_pooled, "pooled")
  results_list[["pooled_het"]] <- extract_tbl(att_pooled_het, "pooled_het")
  results_list[["ai06"]] <- extract_tbl(att_ai06, "ai06")
  if (include_bootstrap) results_list[["bootstrap"]] <- extract_tbl(att_boot, "bootstrap")

  rs <- bind_rows(results_list)

  # ------------------------------------------------------------
  # Add extra rows: pooled_V_E_CI and pooled_het_V_E_CI
  # Replace SE by sqrt(V_E) and build CI using that SE
  # ------------------------------------------------------------
  make_VE_CI_row <- function(df_one_row, new_name) {
    if (is.null(df_one_row) || nrow(df_one_row) == 0) return(NULL)

    # If V_E or att_est is missing, keep NA outputs (still append row)
    ve_se <- ifelse(is.na(df_one_row$V_E), NA_real_, sqrt(df_one_row$V_E))
    att   <- df_one_row$att_est

    df_one_row %>%
      mutate(
        inference_method = new_name,
        SE = ve_se,
        CI_lower = att - 1.96 * ve_se,
        CI_upper = att + 1.96 * ve_se,
        # optional: keep t consistent with this SE
        t = att / ve_se
      )
  }

  pooled_row     <- rs %>% filter(inference_method == "pooled") %>% slice(1)
  pooled_het_row <- rs %>% filter(inference_method == "pooled_het") %>% slice(1)

  pooled_VE_CI_row     <- make_VE_CI_row(pooled_row, "pooled_V_E_CI")
  pooled_het_VE_CI_row <- make_VE_CI_row(pooled_het_row, "pooled_het_V_E_CI")

  rs <- bind_rows(rs, pooled_VE_CI_row, pooled_het_VE_CI_row)


  # True SATT (sample ATT among treated): mean(tau) in sample
  # In gen_df_adv, Y1 = Y0 + tau, so tau = Y1 - Y0 (noise cancels)
  # Here df has Y0, Y1 columns; use treated units Z==1
  true_satt <- df %>%
    filter(Z == 1) %>%
    summarise(att = mean(Y1 - Y0)) %>%
    pull(att)

  rs$SATT <- true_satt

  # Post-match bias: weighted diff in Y0 between treated and controls in matched set
  bias_val <- NA_real_
  full_units <- NULL
  if (!is.null(mtch)) {
    full_units <- tryCatch(result_table(mtch, nonzero_weight_only = TRUE), error = function(e) NULL)
    if (!is.null(full_units) && nrow(full_units) > 0 && all(c("Y0", "weights", "Z") %in% names(full_units))) {
      tmp <- full_units %>%
        group_by(Z) %>%
        summarise(mn = stats::weighted.mean(Y0, w = weights, na.rm = TRUE), .groups = "drop")
      if (nrow(tmp) == 2 && all(c(0, 1) %in% tmp$Z)) {
        bias_val <- tmp$mn[tmp$Z == 1] - tmp$mn[tmp$Z == 0]
      }
    }
  }
  rs$bias <- bias_val

  # annotate metadata
  overlap_tbl <- if (length(overlap_stats) == 0) tibble() else as_tibble(overlap_stats)

  rs <- rs %>%
    bind_cols(overlap_tbl) %>%
    mutate(
      seed = seed_addition,
      gridID = grid_id,
      deg_overlap = overlap_label,
      error_label = error_label,
      prop_nc_unif = prop_nc_unif,
      N = nc + nt,
      nc = nc,
      nt = nt,
      f0_sd = f0_sd,
      ctr_dist = ctr_dist,
      nbins = nbins,
      status = ifelse(is.null(mtch), "Matching Failed", "Success")
    )

  rs
}

# -------------------------------------------------------------------
# sim_master: wrapper called by parallel script
# -------------------------------------------------------------------

sim_master <- function(iteration, N = 600, overlap_label, error_label = "homoskedastic", k_dim = 2, grid_id = 1, ...) {
  # Overlap knob: prop of uniform controls
  prop_nc_unif_values <- c(
    very_high   = 2/3,
    high        = 1/2,
    mid        = 1/3,
    low       = 1/5,
    very_low  = 1/10
  )
  prop_nc_unif <- prop_nc_unif_values[[overlap_label]]
  if (is.null(prop_nc_unif)) stop("Unknown overlap_label: ", overlap_label)

  # fixed sample size per your replication request
  nt <- 100
  nc <- 500

  # seed scheme (stable across grid/labels)
  generate_seed <- function(i, overlap_label, error_label, N) {
    overlap_levels <- c("very_low", "low", "mid", "high", "very_high")
    error_levels <- c("homoskedastic", "covariate_dep", "treatment_dep")
    N_levels <- c(600)

    overlap_index <- match(overlap_label, overlap_levels)
    error_index <- match(error_label, error_levels)
    N_index <- match(N, N_levels)

    if (is.na(overlap_index)) overlap_index <- 1
    if (is.na(error_index)) error_index <- 1
    if (is.na(N_index)) N_index <- 1

    i + 1000 * (overlap_index - 1) + 10000 * (error_index - 1) + 100000 * (N_index - 1)
  }

  seed <- generate_seed(iteration, overlap_label, error_label, N)

  start_time <- Sys.time()
  result <- tryCatch({
    toy_match_infer(
      i = iteration,
      overlap_label = overlap_label,
      error_label = error_label,
      prop_nc_unif = prop_nc_unif,
      nc = nc,
      nt = nt,
      f0_sd = 0.5,
      ctr_dist = 0.5,
      nbins = 6,
      include_bootstrap = TRUE,
      boot_mtd = "wild",
      B = 250,
      seed_addition = seed,
      grid_id = grid_id,
      verbose = FALSE
    )
  }, error = function(e) {
    warning("toy_match_infer failed iter ", iteration, " overlap ", overlap_label, ": ", e$message, call. = FALSE)
    tibble(
      runID = iteration,
      gridID = grid_id,
      deg_overlap = overlap_label,
      error_label = error_label,
      prop_nc_unif = prop_nc_unif,
      N = nc + nt,
      nc = nc,
      nt = nt,
      att_est = NA_real_,
      SE = NA_real_,
      status = "Iteration Failed"
    )
  })
  end_time <- Sys.time()

  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
  result$time_secs <- duration
  result
}
