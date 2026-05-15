# scripts/sims-variance-multi/0_utils.R
#
# Utilities for the 3-factor variance-estimator simulation study.
#
# Factors
# -------
#   1. Overlap (5 levels):        very_low, low, mid, high, very_high
#                                 controlled via prop_nc_unif
#   2. Error structure (2 levels): homo  — σ₀(x) = 0.5  (constant)
#                                  het   — σ₀(x) bimodal Gaussian bumps
#   3. Common variance (2 levels): common    — σ₁ = σ₀  (sigma1_extra = 0)
#                                  no_common — σ₁ = σ₀ + 2 (sigma1_extra = 2)
#   4. k (min number of matches):  1 or 2

# 5 × 2 × 2 x 2 = 40 design cells, 1000 replications.
#
# Three variance estimators
# ------------------------
#   homo       the full homoskedastic
#
#   het        assumes common variance, but allows for heteroskedasticity
#
#   ttmatch    Does treatment-treatment matching to estimate variance on
#.             tx side
#
# Performance measure: coverage of τ_SATT
#   covered = (CI_lower ≤ SATT) & (SATT ≤ CI_upper)


suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(readr)
  library(here)
  library(mvtnorm)
})

# ─────────────────────────────────────────────────────────────────────────────
# Design constants
# ─────────────────────────────────────────────────────────────────────────────

OVERLAP_LABELS <- c("very_low", "low", "mid", "high", "very_high")

PROP_NC_UNIF <- c(
  very_low  = 1/10,
  low       = 1/5,
  mid       = 1/3,
  high      = 1/2,
  very_high = 2/3
)

DESIGN_GRID <- expand_grid(
  overlap_label = OVERLAP_LABELS,
  error_type    = c("homo", "het"),
  sigma1_extra  = c(0, 2),
  k_match = c( 1, 2 )
) %>%
  mutate(
    common_label = if_else(sigma1_extra == 0, "common", "no_common"),
    # factor ordering for plots
    overlap_label = factor(overlap_label, levels = OVERLAP_LABELS),
    error_type    = factor(error_type,    levels = c("homo", "het")),
    common_label  = factor(common_label,  levels = c("common", "no_common"))
  )

# ─────────────────────────────────────────────────────────────────────────────
# σ₀(x): bimodal Gaussian bumps (v2 parameters, A=5, h_sq=0.08, sigma_min=0.2)
# At treated cluster centres (0.25,0.25) and (0.75,0.75): σ₀ ≈ 5.2
# At control cluster centres: σ₀ ≈ 0.2
# ─────────────────────────────────────────────────────────────────────────────

sigma0_bimodal <- function(x1, x2,
                           sigma_min = 0.2, A = 5.0, h_sq = 0.08) {
  d1 <- (x1 - 0.25)^2 + (x2 - 0.25)^2
  d2 <- (x1 - 0.75)^2 + (x2 - 0.75)^2
  sigma_min + A * (exp(-d1 / h_sq) + exp(-d2 / h_sq))
}

# ─────────────────────────────────────────────────────────────────────────────
# Data generation
# ─────────────────────────────────────────────────────────────────────────────

#' Generate one dataset for the multifactor simulation.
#'
#' Uses gen_df_adv_k() with f0_sd = 0 to get the noiseless skeleton
#' (Y0_denoised, Y1_denoised), then adds noise according to error_type
#' and sigma1_extra.
#'
#' @param nc,nt  Control and treated sample sizes.
#' @param prop_nc_unif  Fraction of control units drawn uniformly
#'   (controls degree of overlap).
#' @param ctr_dist  Cluster separation parameter.
#' @param error_type  "homo" (σ₀ = 0.5 constant) or "het" (bimodal σ₀).
#' @param sigma1_extra  Extra treatment-side SD: σ₁ = σ₀ + sigma1_extra.
#'   0 = common variance; 2 = no-common variance.
#' @param seed  Optional seed (set before gen_df_adv_k call).
#' @return Data frame with Y, Y0, Y1, Z, X1, X2, id, sigma0, sigma1.
make_df_multi <- function(
    nc = 500, nt = 100,
    prop_nc_unif,
    ctr_dist     = 0.5,
    error_type   = c("homo", "het"),
    sigma1_extra = 0,
    seed         = NULL) {

  error_type <- match.arg(error_type)
  if (!is.null(seed)) set.seed(seed)

  f0_fun_mat <- function(X) {
    X <- as.matrix(X)
    mvtnorm::dmvnorm(X, mean = c(0.5, 0.5),
                     sigma = matrix(c(1, 0.8, 0.8, 1), 2, 2)) * 20
  }
  tx_fun_mat <- function(X) {
    X <- as.matrix(X)
    3 * X[, 1] + 3 * X[, 2]
  }

  # Generate noiseless skeleton (f0_sd = 0 → Y0_denoised = f0_fun(X))
  df_raw <- gen_df_adv_k(
    nc = nc, nt = nt, k = 2,
    f0_sd         = 0,
    tx_effect_fun = tx_fun_mat,
    f0_fun        = f0_fun_mat,
    ctr_dist      = ctr_dist,
    prop_nc_unif  = prop_nc_unif
  )

  # σ₀ shape depends on error_type factor
  sig0 <- if (error_type == "homo") {
    rep(0.5, nrow(df_raw))
  } else {
    sigma0_bimodal(df_raw$X1, df_raw$X2)
  }

  sig1 <- sig0 + sigma1_extra

  eps0 <- rnorm(nrow(df_raw), 0, sig0)
  # rescale tx errors to match target tx variances
  eps1 <- eps0 * sig1 / sig0

  df_raw %>%
    mutate(
      Y0     = Y0_denoised + eps0,
      Y1     = Y1_denoised + eps1,
      Y      = ifelse(Z, Y1, Y0),
      Z      = as.integer(Z),
      sigma0 = sig0,
      sigma1 = sig1
    )
}



# ─────────────────────────────────────────────────────────────────────────────
# One iteration: match + variance estimators
# ─────────────────────────────────────────────────────────────────────────────

#' Run one simulation iteration for a single design cell.
#'
#' @param i         Replication index (for runID column).
#' @param overlap_label,error_type,sigma1_extra  Design cell identifiers.
#' @param prop_nc_unif  Looked up from PROP_NC_UNIF inside sim_master_multi;
#'   passed directly here.
#' @param seed_addition  Seed (unique per cell × replication).
#' @param K_tt  Treated-to-treated neighbours for ttmatch (default 2).
#' @param verbose  Print error messages from individual estimators.
#' @return Tibble with 4 rows (one per estimator) and all metadata columns.
one_iter <- function(
    i,
    overlap_label,
    error_type    = c("homo", "het"),
    sigma1_extra  = 0,
    prop_nc_unif,
    nc            = 500,
    nt            = 100,
    ctr_dist      = 0.5,
    caliper       = 0.1,
    k_match       = 2,
    K_tt          = 2,
    seed_addition,
    verbose       = FALSE
) {

  error_type <- match.arg(error_type)
  set.seed(seed_addition)

  if ( verbose ) {
    cat( glue::glue("Sim iter #{i}: over:{overlap_label} - error:{error_type} - sigma1:{sigma1_extra} - k={k_match}") )
    cat( "\n" )
  }

  # Make the dataset
  df <- make_df_multi(
    nc = nc, nt = nt,
    prop_nc_unif = prop_nc_unif,
    ctr_dist     = ctr_dist,
    error_type   = error_type,
    sigma1_extra = sigma1_extra,
    seed         = seed_addition
  )

  mtch <-  get_cal_matches(
    data       = df,
    Z ~ X1 + X2,
    rad_method = "adaptive",
    caliper = caliper,
    scaling    = 1,
    k          = k_match,
    warn       = FALSE,
    est_method = "scm"
  )

  true_satt <- df %>%
    filter(Z == 1) %>%
    summarise(att = mean(Y1 - Y0)) %>%
    pull(att)


  ATT <- CSM:::get_att_point_est(mtch,
                                 treatment = "Z",
                                 outcome   = "Y" )

  res_homo <- get_finite_variance(mtch, outcome = "Y",
                                  homoskedastic = TRUE,
                                  use_common_variance = TRUE )

  res_het <- get_finite_variance(mtch, outcome = "Y",
                                 homoskedastic = FALSE,
                                 use_common_variance = TRUE )

  res_tt <- get_finite_variance(mtch, outcome = "Y",
                                homoskedastic = FALSE,
                                use_common_variance = FALSE,
                                K                   = K_tt
  )

  # Now superpopulation ones
  super_homo <- get_total_variance( mtch,
                                    outcome = "Y",
                                    treatment = "Z",
                                    variance_method = "pooled" )

  super_het <- get_total_variance( mtch,
                                   outcome = "Y",
                                   treatment = "Z",
                                   variance_method = "pooled_het" )

  res = bind_rows( homo = res_homo,
                   het = res_het,
                   ttmatch = res_tt,
                   pop_homo = super_homo,
                   pop_het = super_het,
                   .id = "method" )


  res %>%
    mutate( att_est      = ATT,
            SATT         = true_satt,
            runID        = i,
            overlap_label = overlap_label,
            error_type   = error_type,
            sigma1_extra = sigma1_extra,
            common_label = if_else(sigma1_extra == 0, "common", "no_common"),
            nc           = nc,
            nt           = nt,
            k_match      = k_match
    ) %>%
    relocate( runID,
              overlap_label, error_type, common_label, k_match,
              method,
              SATT, att_est )
}



if ( FALSE ) {

  # one_iter() testing code ----
  debugonce(one_iter)
  debugonce()
  one_iter(
    i             = 1,
    overlap_label = "mid",
    error_type    = "homo",
    sigma1_extra  = 2,
    prop_nc_unif  = 0.3,
    nc           = 500,
    nt           = 100,
    seed_addition = 423,
    K_tt         = 2,
    verbose = TRUE
  )


  one_iter(
    i             = 2,
    overlap_label = "low",
    error_type    = "homo",
    sigma1_extra  = 2,
    prop_nc_unif  = 0.3,
    nc           = 500,
    nt           = 100,
    seed_addition = 4243,
    K_tt         = 2,
    verbose = TRUE
  )
}


# ─────────────────────────────────────────────────────────────────────────────
# sim_master_multi: one SLURM job = one iteration × all 20 design cells
# ─────────────────────────────────────────────────────────────────────────────

#' Run replication `iteration` across every row of DESIGN_GRID.
#'
#' Seeds are constructed to be unique per (iteration × cell) while remaining
#' reproducible.  The formula is:
#'   seed = iteration
#'        + 1000  * (overlap_idx - 1)
#'        + 10000 * (error_idx   - 1)   [1=homo, 2=het]
#'        + 50000 * (sigma1_idx  - 1)   [1=common, 2=no_common]
#'
#' @param iteration  Positive integer SLURM array task ID.
#' @param grid  Design grid tibble (default: DESIGN_GRID).
#' @param nc,nt  Sample sizes.
#' @param K_tt  Treated-to-treated neighbours for ttmatch.
#' @return Tibble with 20 × 4 = 80 rows (20 cells × 4 estimators).
sim_master_multi <- function( iteration,
                              grid  = DESIGN_GRID,
                              nc    = 500,
                              nt    = 100,
                              K_tt  = 2,
                              k_match = 2,
                              caliper = 0.1,
                              verbose = FALSE,
                              save_path = NULL,
                              overwrite = FALSE ) {

  if ( !is.null(save_path)  ) {
    out_file <- file.path(save_path, sprintf("iter_%04d.rds", iteration))
    if ( !overwrite && file.exists(out_file) ) {
      return( NA )
    }
  }


  make_seed <- function(iter, overlap_label, error_type, sigma1_extra) {
    overlap_idx <- match(as.character(overlap_label), OVERLAP_LABELS)
    error_idx   <- if (as.character(error_type) == "homo") 1L else 2L
    sigma1_idx  <- if (sigma1_extra == 0) 1L else 2L
    iter +
      1000L  * (overlap_idx - 1L) +
      10000L * (error_idx   - 1L) +
      50000L * (sigma1_idx  - 1L)
  }

  results <- vector("list", nrow(grid))

  for (j in seq_len(nrow(grid))) {

    row  <- grid[j, ]
    seed <- make_seed(iteration,
                      row$overlap_label, row$error_type, row$sigma1_extra)
    prop <- PROP_NC_UNIF[as.character(row$overlap_label)]

    start <- Sys.time()
    results[[j]] <- tryCatch(
      one_iter(
        i             = iteration,
        overlap_label = as.character(row$overlap_label),
        error_type    = as.character(row$error_type),
        sigma1_extra  = row$sigma1_extra,
        prop_nc_unif  = prop,
        caliper = caliper,
        nc            = nc,
        nt            = nt,
        seed_addition = seed,
        K_tt          = K_tt,
        k_match       = row$k_match,
        verbose = verbose
      ),
      error = function(e) {
        warning(sprintf(
          "one_iter failed: iter=%d overlap=%s error=%s sigma1=%.1f — %s",
          iteration, row$overlap_label, row$error_type,
          row$sigma1_extra, e$message
        ), call. = FALSE)
        NA
      }
    )
    results[[j]]$time_secs <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  }

  results <- bind_rows(results)

  if ( !is.null(save_path) ) {
    out_file <- file.path(save_path, sprintf("iter_%04d.rds", iteration))
    cat( "Saving to ", out_file, "  (", nrow(results), " rows, ",
         round(sum(results$time_secs, na.rm = TRUE), 1), "s total)\n", sep = "" )
    saveRDS(results, out_file)
  }

  results
}


if ( FALSE ) {

  sim_master_multi( 1, DESIGN_GRID[1:3,] )
}


# ─────────────────────────────────────────────────────────────────────────────
# Path helpers
# ─────────────────────────────────────────────────────────────────────────────

get_sim_paths <- function(output_name = "sims-variance-multi") {
  base <- here::here("data/outputs", output_name)
  dir.create(file.path(base, "individual"), showWarnings = FALSE, recursive = TRUE)
  list(
    base           = base,
    individual_dir = file.path(base, "individual"),
    combined_rds   = file.path(base, "combined_results.rds"),
    combined_csv   = file.path(base, "combined_results.csv")
  )
}


#' Collect per-iteration .rds files and combine into a single tibble.
#'
#' Reads every file matching \code{iter_NNNN.rds} in \code{individual_dir},
#' binds them into one tibble, and (as a side effect) writes the combined
#' result to \code{combined_results.rds} and \code{combined_results.csv}
#' in the parent output directory.
#'
#' @param output_name  Name of the output subdirectory under
#'   \code{data/outputs/} (default: "sims-variance-multi").  Used to
#'   derive paths via \code{get_sim_paths()}.  Ignored if
#'   \code{individual_dir} is supplied directly.
#' @param individual_dir  Path to the directory containing the
#'   \code{iter_NNNN.rds} files.  If NULL (default), derived from
#'   \code{output_name} via \code{get_sim_paths()}.
#' @param save  If TRUE (default), write the combined tibble to RDS + CSV.
#' @return The combined tibble (invisibly when \code{save = TRUE}).
collect_sim_results <- function(output_name   = "sims-variance-multi",
                                individual_dir = NULL,
                                save           = TRUE) {

  paths <- get_sim_paths(output_name)
  if (is.null(individual_dir)) individual_dir <- paths$individual_dir

  cat("Collecting .rds files from:\n  ", individual_dir, "\n", sep = "")

  files <- list.files(
    path       = individual_dir,
    pattern    = "^iter_\\d+\\.rds$",
    full.names = TRUE
  )

  if (length(files) == 0)
    stop("No iter_NNNN.rds files found in: ", individual_dir)

  cat("Found ", length(files), " file(s).\n", sep = "")

  read_one <- function(fp) {
    out <- tryCatch(readRDS(fp), error = function(e) NULL)
    if (is.null(out)) { warning("Could not read: ", fp); return(NULL) }
    as_tibble(out)
  }

  res_list <- purrr::map(files, read_one)
  failed   <- files[purrr::map_lgl(res_list, is.null)]
  combined <- dplyr::bind_rows(res_list)

  cat("Combined rows : ", nrow(combined), "\n", sep = "")
  cat("Replications  : ", dplyr::n_distinct(combined$runID), "\n", sep = "")
  cat("Design cells  : ", nrow(dplyr::distinct(
    combined, overlap_label, error_type, sigma1_extra, k_match)), "\n", sep = "")

  if (length(failed) > 0) {
    cat("Failed to read ", length(failed), " file(s):\n", sep = "")
    cat(paste("  -", basename(failed), collapse = "\n"), "\n")
  }

  if (save) {
    readr::write_rds(combined, paths$combined_rds)
    readr::write_csv(combined, paths$combined_csv)
    cat("Saved:\n  RDS: ", paths$combined_rds,
        "\n  CSV: ", paths$combined_csv, "\n", sep = "")
    return(invisible(combined))
  }

  combined
}
