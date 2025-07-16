sim_master <- function(iteration, N, overlap_label, error_label, k = 4, grid_id = 1, ...) {
  iteration_results_list <- list()

  # Set overlap level
  prop_nc_unif_values <- c(
    very_low = 2/3,
    low = 1/2,
    mid = 1/3,
    high = 1/5,
    very_high = 1/10
  )
  prop_nc_unif <- prop_nc_unif_values[[overlap_label]]

  # Choose error function
  hetero_fun <- switch(error_label,
                       homoskedastic = function(X, Z) rep(0.5, nrow(X)),
                       covariate_dep = function(X, Z) 0.25 + 0.5 * sqrt(rowSums((X - colMeans(X))^2)),
                       treatment_dep = function(X, Z) ifelse(Z == 1, 0.75, 0.25)
  )

  # Wrapper to generate noise
  noise_fun <- function(X, Z) hetero_fun(X, Z) * rnorm(nrow(X))

  # Set
  nt = round(N / 6)
  nc = N - nt

  start_time <- Sys.time()

  # --- Call one_iteration Function ---
  result <- tryCatch({
    one_iteration(
      i = iteration,
      k = k_dim,
      nc = nc,
      nt = nt,
      prop_nc_unif = prop_nc_unif,
      # true_sigma = true_sigma_base,
      true_sigma = NULL,  # handled by noise_fun
      noise_fun = noise_fun,  # needs to be added to one_iteration
      # supporting parameters
      toy_ctr_dist = toy_ctr_dist_base,
      scaling = scaling_base,
      include_bootstrap = include_bootstrap_global,
      boot_mtd = boot_mtd_global,
      B = B_global,
      seed_addition = iteration,
      verbose = FALSE
    )
  }, error = function(e) {
    warning("one_iteration failed for iteration ", iteration, ", overlap ", overlap_label, ": ", e$message, call. = FALSE)
    log_path <- file.path("/homes2/xmeng/scmatch2/scripts/inference-scripts/logs", sprintf("error_iter%d_grid%d.txt", iteration, grid_id))
    dir.create("logs", showWarnings = FALSE, recursive = TRUE)

    writeLines(c(
      sprintf("Timestamp: %s", Sys.time()),
      sprintf("Iteration: %d", iteration),
      sprintf("Grid ID: %d", grid_id),
      sprintf("N: %d, Overlap: %s, Error: %s", N, overlap_label, error_label),
      sprintf("R library path is: %s", .libPaths()),
      "Error message:",
      conditionMessage(e),
      "Traceback:",
      paste(capture.output(traceback()), collapse = "\n")
    ), con = log_path)


    return(tibble(
      runID = iteration, # Assuming one_iteration would add this, but add here for consistency
      k = k_dim,
      # deg_overlap = overlap_label, # no need to add here because can add sepeareatly below
      # prop_nc_unif = prop_nc_unif,
      N = N,
      nc = nc,
      nt = nt,
      error_label = error_label,
      toy_ctr_dist = toy_ctr_dist_base, # Add toy_ctr_dist here
      scaling = scaling_base,         # Add scaling here
      # true_sigma = true_sigma_base,   # Add true_sigma here
      SATT = NA_real_,            # Placeholder for Sample ATT
      att_est = NA_real_,
      SE = NA_real_,
      CI_lower = NA_real_,
      CI_upper = NA_real_,
      bias = NA_real_,            # Placeholder for bias
      N_T = NA_integer_,
      ESS_C = NA_real_,
      V_E = NA_real_,
      sigma_hat = NA_real_,
      status = "Iteration Failed" # Add a status column
      # Add other columns returned by one_iteration as NA if needed
    ))
  }
  )

  # --- Record Time and Store Results ---
  end_time <- Sys.time()
  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
  # cat("      Overlap", overlap_label, "took", round(duration, 2), "seconds.\n")


  # -- Annotate metadata
  result %>%
    dplyr::mutate(
      runID = iteration,
      gridID = grid_id,
      deg_overlap = overlap_label,
      error_label = error_label,
      prop_nc_unif = prop_nc_unif,
      N = N,
      nc = nc,
      nt = nt,
      time_secs = duration,
      status = ifelse("status" %in% names(.), status, "Success")
    )

}

#' Run one iteration of the simulation
#'
#' @param i Iteration number (integer).
#' @param k Dimensionality of covariates (integer, default 2).
#' @param toy_ctr_dist Control cluster distance for gen_one_toy (numeric).
#' @param scaling Scaling parameter for get_cal_matches (numeric).
#' @param true_sigma SD of the noise added to potential outcomes (numeric).
#' @param noise_fun Noise function specification (function of (covar, treatment variable))
#' @param prop_nc_unif Proportion of uniform controls for gen_one_toy (numeric).
#' @param nc Number of control units (integer).
#' @param nt Number of treated units (integer, passed to gen_one_toy).
#' @param save_match Save the resulting match object? (logical, default FALSE).
#' @param path_to_save_match Directory path to save the match object if save_match is TRUE (character). Filename will be generated based on iteration number.
#' @param include_bootstrap Include bootstrap inference results (logical, default TRUE).
#' @param boot_mtd Bootstrap method when include_bootstrap = TRUE. Options: "wild", "Bayesian", "naive" (character, default "wild").
#' @param B Number of bootstrap samples (integer, default 250).
#' @param seed_addition Additional seed for bootstrap reproducibility (integer, default 11).
#' @param verbose Print progress message (logical).
#' @param R Total number of iterations (for verbose printing).
#'
#' @return A tibble with results for one iteration, including both pooled and bootstrap results (if requested), plus overlap statistics.
#' @importFrom stats rnorm weighted.mean reformulate quantile
#' @importFrom dplyr %>% rename mutate select filter summarize pull group_by last first any_of n tibble bind_rows

one_iteration <- function( i,
                           k = 2,
                           toy_ctr_dist = 0.5,
                           scaling = 8,
                           true_sigma = 0.5,
                           noise_fun = NULL,
                           prop_nc_unif = 1/3,
                           nc = 500,
                           nt = 100,
                           save_match = FALSE,
                           path_to_save_match = "",
                           include_bootstrap = TRUE,
                           boot_mtd = "wild",
                           B = 250,
                           seed_addition = 11,
                           verbose = TRUE,
                           R = NA ) {

  # # Verbose output
  # if( verbose ) {
  #   if (!is.na(R) && R > 0 && (i %% 10 == 0 || (i == R)) ) {
  #     cat("iteration", i, " (k=", k, ")\n", sep="")
  #   } else if (is.na(R) && i %% 10 == 0) {
  #     cat("iteration", i, " (k=", k, ")\n", sep="")
  #   }
  # }

  ### Make the data using k dimensions
  df_dgp <-
    gen_one_toy( k = k,
                 nc = nc,
                 nt = nt,
                 ctr_dist=toy_ctr_dist,
                 prop_nc_unif=prop_nc_unif,
                 f0_sd = 0 # (optional) set noise level to zero
                 )

  ## Add in the desired noise level for this iteration
  df_dgp_i <- df_dgp %>%
    mutate(noise = if (!is.null(noise_fun)) noise_fun(select(., starts_with("X")), Z) else rnorm(n(), 0, true_sigma)) %>%
    # mutate(noise = rnorm(n(), mean=0, sd=true_sigma)) %>%
    mutate(Y0 = Y0_denoised + noise,
           Y1 = Y1_denoised + noise,
           Y = ifelse(Z, Y1, Y0))



  ### Perform matching
  mtch <- tryCatch(
    get_cal_matches(
      df_dgp_i,
      metric = "maximum",
      scaling = scaling,
      caliper = 1,
      rad_method = "adaptive-5nn",
      est_method = "scm"
    ), error = function(e) {
      warning("Matching failed in iteration ", i, ": ", e$message, call. = FALSE)
      return(NULL)
    }
  )

  # Save the match object if requested
  if (save_match && !is.null(mtch)) {
    if (nzchar(path_to_save_match) && dir.exists(path_to_save_match)) {
      file_name <- paste0(
        "match_iter_", i,
        "_k_", k,
        "_dist_", toy_ctr_dist,
        "_scale_", scaling,
        "_sigma_", true_sigma,
        "_prop_", round(prop_nc_unif, 4),
        "_nc_", nc,
        "_nt_", nt,
        ".rds"
      )
      full_path <- file.path(path_to_save_match, file_name)

      tryCatch({
        saveRDS(mtch, file = full_path)
        if (verbose) {
          cat("   -> Saved match object for iteration", i, "to:", full_path, "\n")
        }
      }, error = function(e) {
        warning("Failed to save match object for iteration ", i, " to ", full_path, ": ", e$message, call. = FALSE)
      })
    } else if (nzchar(path_to_save_match)) {
      warning("Directory '", path_to_save_match, "' does not exist. Cannot save match object for iteration ", i, ".", call. = FALSE)
    } else {
      warning("'path_to_save_match' is empty. Cannot save match object for iteration ", i, ".", call. = FALSE)
    }
  } else if (save_match && is.null(mtch)) {
    if (verbose) {
      cat("   -> Matching failed for iteration", i, ", skipping save.\n")
    }
  }

  ### Calculate overlap statistics
  overlap_stats <- list(
    avg_shared_controls = NA,
    p75_shared_controls = NA,
    avg_shared_treated = NA,
    p75_shared_treated = NA
  )

  if (!is.null(mtch)) {
    overlap_stats <- tryCatch({
      calculate_overlap_statistics(mtch)
    }, error = function(e) {
      warning("Overlap statistics calculation failed in iteration ", i, ": ", e$message, call. = FALSE)
      list(
        avg_shared_controls = NA,
        p75_shared_controls = NA,
        avg_shared_treated = NA,
        p75_shared_treated = NA
      )
    })
  }

  ### Perform inference - both pooled and bootstrap (if requested)
  results_list <- list()

  # Always do pooled inference
  ATT_estimate_pooled <- NULL
  if (!is.null(mtch)) {
    ATT_estimate_pooled <- tryCatch(
      get_ATT_estimate( mtch, variance_method = "pooled" ),
      error = function(e) {
        warning("Pooled ATT estimation failed in iteration ", i, ": ", e$message, call. = FALSE)
        return(NULL)
      }
    )
  }

  # Add bootstrap inference if requested
  ATT_estimate_bootstrap <- NULL
  if (include_bootstrap && !is.null(mtch)) {
    ATT_estimate_bootstrap <- tryCatch(
      get_ATT_estimate( mtch,
                        variance_method = "bootstrap",
                        boot_mtd = boot_mtd,
                        B = B,
                        seed_addition = seed_addition ),
      error = function(e) {
        warning("Bootstrap ATT estimation failed in iteration ", i, ": ", e$message, call. = FALSE)
        return(NULL)
      }
    )
  }

  # Helper function to extract results safely
  extract_results <- function(ATT_estimate, method_name) {
    att_est_val <- ATT_estimate$ATT %||% NA
    se_val <- ATT_estimate$SE %||% NA
    N_T_val <- ATT_estimate$N_T %||% NA
    ESS_C_val <- ATT_estimate$ESS_C %||% NA
    V_E_val <- ATT_estimate$V_E %||% NA
    V_P_val <- ATT_estimate$V_P %||% NA
    sigma_hat_val <- ATT_estimate$sigma_hat %||% NA

    tibble(
      runID = i,
      k = k,
      inference_method = method_name,
      boot_mtd = ifelse(method_name == "bootstrap", boot_mtd, NA_character_),
      B = ifelse(method_name == "bootstrap", B, NA_integer_),
      att_est = att_est_val,
      SE = se_val,
      CI_lower = att_est_val - 1.96 * se_val,
      CI_upper = att_est_val + 1.96 * se_val,
      N_T = N_T_val,
      ESS_C = ESS_C_val,
      V_E = V_E_val,
      V_P = V_P_val,
      sigma_hat = sigma_hat_val,
      # Add overlap statistics
      avg_shared_controls = overlap_stats$avg_shared_controls,
      p75_shared_controls = overlap_stats$p75_shared_controls,
      avg_shared_treated = overlap_stats$avg_shared_treated,
      p75_shared_treated = overlap_stats$p75_shared_treated
    )
  }

  # Extract pooled results
  if (!is.null(ATT_estimate_pooled)) {
    results_list[["pooled"]] <- extract_results(ATT_estimate_pooled, "pooled")
  } else {
    # Create empty row for failed pooled estimation
    results_list[["pooled"]] <- tibble(
      runID = i, k = k, inference_method = "pooled",
      boot_mtd = NA_character_, B = NA_integer_,
      att_est = NA, SE = NA, CI_lower = NA, CI_upper = NA,
      N_T = NA, ESS_C = NA, V_E = NA, V_P = NA, sigma_hat = NA,
      avg_shared_controls = NA,
      p75_shared_controls = NA,
      avg_shared_treated = NA,
      p75_shared_treated = NA
    )
  }

  # Extract bootstrap results if requested
  if (include_bootstrap) {
    if (!is.null(ATT_estimate_bootstrap)) {
      results_list[["bootstrap"]] <- extract_results(ATT_estimate_bootstrap, "bootstrap")
    } else {
      # Create empty row for failed bootstrap estimation
      results_list[["bootstrap"]] <- tibble(
        runID = i, k = k, inference_method = "bootstrap",
        boot_mtd = boot_mtd, B = B,
        att_est = NA, SE = NA, CI_lower = NA, CI_upper = NA,
        N_T = NA, ESS_C = NA, V_E = NA, V_P = NA, sigma_hat = NA,
        avg_shared_controls = NA,
        p75_shared_controls = NA,
        avg_shared_treated = NA,
        p75_shared_treated = NA
      )
    }
  }

  # Combine results
  rs <- bind_rows(results_list)

  # Get true ATT (based on denoised potential outcomes for treated)
  true_att <- df_dgp %>%
    filter(Z == TRUE) %>%
    summarize(att = mean(Y1_denoised - Y0_denoised)) %>%
    pull(att)

  rs$SATT <- true_att

  # Calculate true post-match bias (same for both methods)
  bias_val <- NA
  if (!is.null(mtch)) {
    full_units <- tryCatch({
      full_unit_table(mtch, nonzero_weight_only = TRUE)
    }, error = function(e) {
      warning("Failed to get full unit table in iteration ", i, ": ", e$message, call. = FALSE)
      return(NULL)
    })

    if (!is.null(full_units) && nrow(full_units) > 0 &&
        "Y0_denoised" %in% names(full_units) && "weights" %in% names(full_units)) {
      bias_calc <- full_units %>%
        group_by(Z) %>%
        summarize(mn = stats::weighted.mean(Y0_denoised, w = weights, na.rm=TRUE),
                  .groups = 'drop')

      if(nrow(bias_calc) == 2 && all(c(0, 1) %in% bias_calc$Z)) {
        mean_y0_t_matched <- bias_calc$mn[bias_calc$Z == 1]
        mean_y0_c_matched <- bias_calc$mn[bias_calc$Z == 0]
        bias_val <- mean_y0_t_matched - mean_y0_c_matched
      } else {
        warning("Bias calculation skipped in iteration ", i, ": Both Z groups not present or valid in weighted data.", call. = FALSE)
      }
    } else if (!is.null(full_units)) {
      if (!"Y0_denoised" %in% names(full_units)) warning("Bias calculation skipped in iteration ", i, ": Y0_denoised not found in full_units.", call. = FALSE)
      if (!"weights" %in% names(full_units)) warning("Bias calculation skipped in iteration ", i, ": 'weights' not found in full_units.", call. = FALSE)
    }
  }

  rs$bias <- bias_val

  return(rs)
}
#' Get Matched Matrix
#'
#' Converts a full matched table into a matrix where each row corresponds to a treated unit
#' and each column contains the IDs of matched control units. Missing values (if any) are filled with `NA`.
#' Fixed version that properly handles control reuse across subclasses.
#'
#' @param full_matched_table A data frame containing the matched data. Must include the columns:
#' \itemize{
#'   \item{\code{id}: Unique identifiers for each unit.}
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
#'   id = 1:6,
#'   Z = c(1, 0, 1, 0, 0, 1),
#'   subclass = c(1, 1, 2, 2, 2, 3)
#' )
#' get_matched_matrix(full_matched_table)
#'
#' @export
get_matched_matrix <- function(full_matched_table) {
  # Get unique treated units
  treated_units <- unique(full_matched_table$id[full_matched_table$Z == 1])

  matched_control_ids <- list()

  # For each treated unit, get ALL controls it's matched to (across all subclasses)
  for (treated_id in treated_units) {
    # Find all subclasses this treated unit appears in
    treated_subclasses <- unique(full_matched_table$subclass[full_matched_table$id == treated_id & full_matched_table$Z == 1])

    # Get all controls from those subclasses
    all_controls <- c()
    for (subclass in treated_subclasses) {
      controls_in_subclass <- full_matched_table$id[full_matched_table$subclass == subclass & full_matched_table$Z == 0]
      all_controls <- c(all_controls, controls_in_subclass)
    }

    # Remove duplicates - this is the key fix!
    # Convert to numeric if possible for proper handling
    if (all(grepl("^[0-9]+$", all_controls))) {
      unique_controls <- unique(as.numeric(all_controls))
    } else {
      unique_controls <- unique(all_controls)
    }

    matched_control_ids[[as.character(treated_id)]] <- unique_controls
  }

  # Create matrix
  max_controls <- max(sapply(matched_control_ids, length))
  matched_matrix <- do.call(rbind, lapply(matched_control_ids, function(ids) {
    c(ids, rep(NA, max_controls - length(ids)))
  }))

  # Ensure proper row names (treated unit IDs)
  rownames(matched_matrix) <- names(matched_control_ids)

  return(matched_matrix)
}

#' Compute Shared Neighbors
#'
#' Computes the number of shared control neighbors between each pair of treated units
#' based on a given matched matrix. Fixed version that properly handles numeric IDs and NA values.
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

  # Initialize with zeros (not FALSE as in original)
  shared_neighbors <- matrix(0, nrow = N1, ncol = N1)
  rownames(shared_neighbors) <- rownames(matched_matrix)
  colnames(shared_neighbors) <- rownames(matched_matrix)

  for (i in 1:(N1 - 1)) {
    for (j in (i + 1):N1) {
      # Get non-NA controls for each treated unit
      controls_i <- matched_matrix[i, ]
      controls_j <- matched_matrix[j, ]

      # Remove NAs before computing intersection
      controls_i <- controls_i[!is.na(controls_i)]
      controls_j <- controls_j[!is.na(controls_j)]

      # Find intersection
      shared <- intersect(controls_i, controls_j)
      shared_count <- length(shared)

      shared_neighbors[i, j] <- shared_neighbors[j, i] <- shared_count
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
#' @param shared_neighbors A numeric matrix of dimensions \eqn{N1 \times N1}, where each entry represents
#' the number of shared neighbors between treated units i and j.
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
#' shared_neighbors <- matrix(c(0, 1, 2, 1, 0, 1, 2, 1, 0), nrow = 3, byrow = TRUE)
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

# Helper function to calculate overlap statistics following the 3-step process
calculate_overlap_statistics <- function(mtch) {
  # Get the full matched table from the match object
  full_matched_table <- tryCatch({
    full_unit_table(mtch, nonzero_weight_only = TRUE)
  }, error = function(e) {
    warning("Failed to get full unit table for overlap statistics: ", e$message, call. = FALSE)
    return(NULL)
  })

  if (is.null(full_matched_table) || nrow(full_matched_table) == 0) {
    return(list(
      avg_shared_controls = NA,
      p75_shared_controls = NA,
      avg_shared_treated = NA,
      p75_shared_treated = NA
    ))
  }

  # Step 1: Get matched matrix (now uses fixed version)
  matched_matrix <- tryCatch({
    get_matched_matrix(full_matched_table)
  }, error = function(e) {
    warning("Failed to get matched matrix for overlap statistics: ", e$message, call. = FALSE)
    return(NULL)
  })

  if (is.null(matched_matrix)) {
    return(list(
      avg_shared_controls = NA,
      p75_shared_controls = NA,
      avg_shared_treated = NA,
      p75_shared_treated = NA
    ))
  }

  # Step 2: Compute shared neighbors (now uses fixed version)
  shared_neighbors <- tryCatch({
    compute_shared_neighbors(matched_matrix)
  }, error = function(e) {
    warning("Failed to compute shared neighbors for overlap statistics: ", e$message, call. = FALSE)
    return(NULL)
  })

  if (is.null(shared_neighbors)) {
    return(list(
      avg_shared_controls = NA,
      p75_shared_controls = NA,
      avg_shared_treated = NA,
      p75_shared_treated = NA
    ))
  }

  # Step 3: Compute overlap statistics (unchanged)
  overlap_stats <- tryCatch({
    compute_overlap_statistics(shared_neighbors)
  }, error = function(e) {
    warning("Failed to compute overlap statistics: ", e$message, call. = FALSE)
    return(list(
      avg_shared_controls = NA,
      p75_shared_controls = NA,
      avg_shared_treated = NA,
      p75_shared_treated = NA
    ))
  })

  return(overlap_stats)
}

# Utility function for safe requireNamespace/library loading
load_libs <- function(libs) {
  for (lib in libs) {
    if (!requireNamespace(lib, quietly = TRUE)) {
      warning(paste("Package", lib, "not found. Please install it."), call. = FALSE)
      return(FALSE)
    }
    library(lib, character.only = TRUE)
  }
  return(TRUE)
}

# `%||%` operator for safe default values
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (is.atomic(x) && all(is.na(x)))) y else x
}



kurtosis <- function(x){
  S_T = sd(x)
  kurt = mean( (x - mean(x))^4 ) / S_T^4
  return(kurt)
}

MCvar_SE <- function(x){
  S_T = sd(x); R <- length(x); k_T <- kurtosis(x)
  return(S_T^2 * sqrt( (k_T-1)/R ))
}

MCSE_SE <- function(x){
  return(sqrt(MCvar_SE(x)))
}

MCSE_bias <- function(x){
  S_T = sd(x); R <- length(x)
  return( sqrt(S_T^2 /R ))
}




#' Save results to file, possibly appending them if file exists.
save_res_to_csv<- function(curr_res,
                           FNAME) {

  if (file.exists(FNAME)) {
    write_csv(curr_res, FNAME, append=TRUE)
  } else {
    write_csv(curr_res, FNAME)
  }

} # save_res_to_csv


get_sim_paths <- function() {
  output_dir_base <- here::here("data/outputs/inf_toy")
  list(
    output_dir_base = output_dir_base,
    individual_dir = file.path(output_dir_base, "individual"),
    summary_csv = file.path(output_dir_base, "summary_table.csv"),
    combined_csv = file.path(output_dir_base, "combined_results.csv")
  )
}


# --- Collection Script for 4D Bias Inference Simulation Results ---
# This script reads all result files from simulation iterations and combines them
# into a single stacked CSV file.
collect_results_to_csv <- function(){
  library(tidyverse)
  library(here)

  # --- Configuration ---
  # # Base directory (same as in the simulation script)
  # output_dir_base <- here::here("data/outputs/4d_bias_inference_Apr2025")
  # individual_dir <- file.path(output_dir_base, "estimation")
  # combined_csv_file <- file.path(output_dir_base, "combined_results.csv")
  paths <- get_sim_paths()
  output_dir_base <- paths$output_dir_base
  individual_dir <- paths$individual_dir
  combined_csv_file <- paths$combined_csv


  # --- Get list of all result files ---
  cat("Looking for result files in:", individual_dir, "\n")

  result_files <- list.files(
    path = individual_dir,
    pattern = "^results_iter_\\d+\\.rds$",
    full.names = TRUE
  )

  num_files <- length(result_files)

  if (num_files == 0) {
    stop("No result files found in the directory. Check path and file pattern.", call. = FALSE)
  } else {
    cat("Found", num_files, "result files.\n")
  }

  # --- Read and combine all results ---
  cat("Reading and combining result files...\n")

  # Progress tracking
  total_files <- length(result_files)
  progress_step <- max(1, floor(total_files / 10))  # Report progress at 10% increments

  # Read all files and combine into a single dataframe
  combined_results <- tibble()
  failed_files <- character()

  for (i in seq_along(result_files)) {
    file_path <- result_files[i]

    # Report progress
    if (i %% progress_step == 0 || i == total_files) {
      percent_done <- round(i / total_files * 100)
      cat(sprintf("  Progress: %d%% (%d/%d files)\n", percent_done, i, total_files))
    }

    # Try to read the file
    result <- tryCatch({
      readRDS(file_path)
    }, error = function(e) {
      file_name <- basename(file_path)
      warning("Failed to read file: ", file_name, " - ", e$message, call. = FALSE)
      failed_files <- c(failed_files, file_name)
      return(NULL)
    })

    # If successful, add to combined results
    if (!is.null(result)) {
      combined_results <- bind_rows(combined_results, result)
    }
  }

  # --- Report on the combined dataset ---
  if (nrow(combined_results) == 0) {
    stop("No valid results were found. Check if files are properly formatted.", call. = FALSE)
  } else {
    cat("\nSuccessfully combined", nrow(combined_results), "rows from",
        total_files - length(failed_files), "files.\n")

    if (length(failed_files) > 0) {
      cat("Failed to read", length(failed_files), "files.\n")
    }

    # Summary of iterations and overlap degrees
    iterations <- length(unique(combined_results$runID))
    overlap_degrees <- unique(combined_results$deg_overlap)

    cat("Dataset contains", iterations, "iterations across",
        length(overlap_degrees), "overlap degrees:",
        paste(overlap_degrees, collapse=", "), "\n")
  }

  # --- Save combined results to CSV ---
  cat("Saving combined results to:", combined_csv_file, "\n")

  tryCatch({
    write_csv(combined_results, combined_csv_file)
    cat("Combined results successfully saved.\n")
  }, error = function(e) {
    stop("Failed to save the combined CSV file: ", e$message, call. = FALSE)
  })

  # # --- Create a summary statistics file ---
  # summary_file <- file.path(output_dir_base, "results_summary.csv")
  # cat("Generating summary statistics...\n")
  #
  # # Calculate summary statistics by overlap degree
  # summary_stats <- combined_results %>%
  #   mutate(ATT = 1.5 * k) %>%
  #   group_by(deg_overlap) %>%
  #   summarize(
  #     n_iterations = n_distinct(runID),
  #     mean_bias = mean(bias, na.rm = TRUE),
  #     sd_bias = sd(bias, na.rm = TRUE),
  #     mean_SATT = mean(SATT, na.rm = TRUE),
  #     mean_att_est = mean(att_est, na.rm = TRUE),
  #     mean_SE = mean(SE, na.rm = TRUE),
  #     coverage_rate = mean(CI_lower <= ATT & ATT <= CI_upper, na.rm = TRUE),
  #     coverage_rate_debiased = mean(CI_lower - bias <= ATT & ATT <= CI_upper-bias, na.rm = TRUE),
  #     mean_ESS_C = mean(ESS_C, na.rm = TRUE),
  #     mean_V_E = mean(V_E,na_rm=T),
  #     mean_sigma_hat = mean(sigma_hat, na.rm = TRUE),
  #     median_time_secs = median(time_secs, na.rm = TRUE),
  #     success_rate = mean(status == "Success", na.rm = TRUE)
  #   )
  #
  # # Save summary statistics
  # tryCatch({
  #   write_csv(summary_stats, summary_file)
  #   cat("Summary statistics saved to:", summary_file, "\n")
  # }, error = function(e) {
  #   warning("Failed to save summary statistics: ", e$message, call. = FALSE)
  # })

  cat("\nCollection process complete.\n")
}
