library(MASS) # For mvrnorm
library(dplyr) # For data manipulation
library(here)
library(CSM)

generate_dgp <- function(N1, N0, panel = "A") {
  X <- c(rnorm(N1, mean = 0, sd = 1),
         rnorm(N0, mean = 0, sd = 1))

  Z <- c(rep(1, N1), rep(0, N0))

  if (panel == "A") {
    mu_1 <- function(x) x
    epsilon_1 <- function() rnorm(1, mean = 0, sd = 1)
  } else if (panel == "B") {
    mu_1 <- function(x) (qchisq(pnorm(x), df = 8) - 8) / sqrt(16)
    epsilon_1 <- function() rnorm(1, mean = 0, sd = 1)
  } else if (panel == "C") {
    mu_1 <- function(x) (qchisq(pnorm(x), df = 1) - 1) / sqrt(2)
    epsilon_1 <- function() rnorm(1, mean = 0, sd = 1)
  } else if (panel == "D") {
    mu_1 <- function(x) (qchisq(pnorm(x), df = 1) - 1) / sqrt(2)
    epsilon_1 <- function() (rchisq(1, df = 1) - 1) / sqrt(2)
  } else if (panel == "E") {
    mu_1 <- function(x) (qchisq(pnorm(x), df = 1) - 1) / sqrt(2)
    epsilon_1 <- function() 2 * (rchisq(1, df = 1) - 1) / sqrt(2)
  }else if (panel == "F") {
    # Panel F: Distribution of tau(x)
    mu_1 <- function(x) {
      10 * sin(pi * x[1] * x[2]) +
        20 * (x[3] - 0.5)^2 +
        10 * x[4] + 5 * x[5] +
        x[3] * cos(pi * x[1] * x[2])
    }
    epsilon_1 <- function() rnorm(1, mean = 0, sd = 1)
    generate_X <- function(n) matrix(rnorm(n * 10), nrow = n, ncol = 10)
  } else {
    stop("Invalid panel specified. Choose from 'A', 'B', 'C', 'D', 'E', or 'F'.")
  }

  epsilon_0_values <- rnorm(N1 + N0, mean = 0, sd = 1)
  if (panel == "F") {
    X <- generate_X(N1 + N0)
    tau_x <- X[, 3] * cos(pi * X[, 1] * X[, 2])
    mu_1_values <- apply(X, 1, mu_1)
    mu_0_values <- mu_1_values - tau_x
    epsilon_1_values <- rnorm(N1 + N0, mean = 0, sd = 1)
  } else {
    epsilon_1_values <- sapply(1:(N1 + N0), function(i) epsilon_1())
    mu_1_values <- mu_1(X)
    mu_0_values <- rep(0, N1 + N0)
  }

  # Generate outcomes
  Y1 <- mu_1_values + epsilon_1_values
  Y0 <- mu_0_values + epsilon_0_values
  Y <- ifelse(Z == 1, Y1, Y0)
  data <- data.frame(
    ID = 1:(N1 + N0),
    X = X,
    Z = Z,
    Y0 = Y0,
    Y1 = Y1,
    Y = Y,
    noise_0 = epsilon_0_values,
    noise_1 = epsilon_1_values
  )

  return(data)
}


## Next: implement KNN method in get_cal_matches function

generate_full_matched_table <- function(dat,
         M,
         file_name = "one-full-table.rds"){
  # devtools::load_all()
  mtch <- get_cal_matches(
    dat,
    # covs = "X",
    covs = starts_with("X"),
    treatment = "Z",
    scaling = 1,
    metric = "euclidean",
    rad_method = "knn",
    k = M
  )

  full_matched_table <- full_unit_table(mtch,
                                nonzero_weight_only = F )
  saveRDS(full_matched_table,
          file = here("scripts/new-inference/data/",file_name))

}



# ## The below function took into the original DGP code
# ## Original DGP format: (unmatched )
# ##    treated = data.frame(ID = 1:N1, X = X_treated, Y = Y_treated),
# #     control = data.frame(ID = 1:N0, X = X_control, Y = Y_control)
# ## Current data format: name is full_table
# ##  it is a matched dataset where each row is just one data entry
# ##  variable Z is used to tell treated or control
# ## variable "subclass" is the identifier of the treated and the matched controls
# ## so the matched controls will have the same subclass value of the treated
# ## i want the function to take into full_table and output the same
# ## matched_pairs variable
# match_controls <- function(treated, control, M) {
#   matched_pairs <- lapply(1:nrow(treated), function(i) {
#     distances <- abs(control$X - treated$X[i])
#     matched_indices <- order(distances)[1:M]
#     list(
#       controls = control$Y[matched_indices],
#       control_indices = matched_indices
#     )
#   })
#   return(matched_pairs)
# }
#

get_matched_control_ids <- function(full_matched_table) {

  subclasses <- unique(full_matched_table$subclass)

  matched_control_ids <- list()

  for (subclass in subclasses) {
    data_in_subclass <- full_matched_table[full_matched_table$subclass == subclass, ]

    treated_units <- data_in_subclass[data_in_subclass$Z == 1, ]
    control_units <- data_in_subclass[data_in_subclass$Z == 0, ]

    for (treated_id in treated_units$ID) {
      matched_control_ids[[as.character(treated_id)]] <- control_units$ID
    }
  }

  max_controls <- max(sapply(matched_control_ids, length))
  matched_matrix <- do.call(rbind, lapply(matched_control_ids, function(ids) {
    c(ids, rep(NA, max_controls - length(ids)))
  }))

  rownames(matched_matrix) <- names(matched_control_ids)

  return(matched_matrix)
}

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


compute_tau_i <- function(full_matched_table,
                          tau_0 = 0){
  N1 <- length(unique(full_matched_table$subclass))
  tau_i <- full_matched_table %>%
    group_by(Z, subclass) %>%
    summarize(avg_Y = sum(Y*weights) / sum(weights)) %>%
    group_by(subclass) %>%
    summarize(est = last(avg_Y) - first(avg_Y)) %>%
    pull(est) - tau_0
  return(tau_i)
}

compute_t_stat <- function(S, sd_method = "usual",
                           signs = NULL,
                           full_matched_table = NULL){
  if (sd_method == "usual"){
    sd_S = sd(S)
  }else if (sd_method == "pooled") {
    weighted_var <- get_pooled_variance(
      matches_table = full_matched_table,
      outcome = "Y",
      treatment = "Z",
      var_weight_type = "uniform"
    )
    sigma_hat <- sqrt(weighted_var)

    unique_subclasses <- full_matched_table %>%
      filter(Z == TRUE) %>%
      distinct(subclass) %>%
      arrange(subclass) %>% # Ensure a stable ordering
      mutate(sign = signs)

    full_matched_table_signed <- full_matched_table %>%
      left_join(unique_subclasses, by = "subclass") %>%
      mutate(weights = weights * sign) %>%
      select(-sign)

    Ns_signed <- calc_N_T_N_C(full_matched_table_signed)

    sd_S <- get_plug_in_SE(
      N_T = Ns_signed$N_T,
      ESS_C = Ns_signed$N_C_tilde,
      sigma_hat = sigma_hat
    )
  }else{
    cat("sd_method in compute_t_stat should be either `usual` or `pooled`" )
  }
  t_stat <- abs(mean(S)) / sd_S
  return( t_stat )
}

ferman_sign_change_null_dist <-
  function(full_matched_table,
           tau_0 = 0,
           max_permutations = 1000,
           overlap_adj = T,
           sd_method = "usual"){
    N1 <- length(unique(full_matched_table$subclass))

    matched_matrix <- get_matched_control_ids(full_matched_table = full_matched_table)
    shared_neighbors <- compute_shared_neighbors(matched_matrix)
    shared_neighbors_binary <- shared_neighbors > 0


    sample_signs <- function(N1, shared_neighbors_binary, overlap_adj){
      signs <- numeric(N1)
      if (overlap_adj){
        signs[1] <- sample(c(-1, 1), 1)
        for (i in 2:N1) {
          connected <- which(shared_neighbors_binary[i, 1:(i - 1)])

          if (length(connected) > 0) {
            signs[i] <- signs[connected[1]]
          } else {
            signs[i] <- sample(c(-1, 1), 1)
          }
        }
      }else{
        signs <- sample(c(-1, 1), N1, replace=T)
      }
      return(signs)
    }

    T_null <- replicate(max_permutations, {
      signs <- sample_signs(
        N1=N1,
        shared_neighbors_binary = shared_neighbors_binary,
        overlap_adj = overlap_adj)
      tau_i <- compute_tau_i(full_matched_table,
                             tau_0)
      gS <- signs * tau_i
      compute_t_stat(gS,
                     sd_method,
                     signs,
                     full_matched_table)
    })
    return(T_null)
  }

ferman_sign_change_test <-
  function(full_matched_table,
           tau_0 = 0,
           alpha = 0.05,
           max_permutations = 1000,
           overlap_adj = T,
           sd_method = "usual") {
    tau_i <- compute_tau_i(full_matched_table,
                           tau_0)
    T_obs <- compute_t_stat(tau_i)
    T_null <- ferman_sign_change_null_dist(
      full_matched_table = full_matched_table,
      tau_0 = tau_0,
      max_permutations = max_permutations,
      overlap_adj = overlap_adj,
      sd_method = sd_method)
    critical_value <-
      quantile(T_null, probs = 1 - alpha)
    reject <-
      as.numeric(T_obs > critical_value)
    reject
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


# compute_overlap_statistics_old <- function(matched_pairs) {
#   N1 <- length(matched_pairs)
#
#   shared_controls <- numeric(N1) # Initialize a vector to store shared control counts
#   shared_treated <- numeric(N1)  # Initialize a vector to store shared treated counts
#
#   for (i in 1:N1) {
#     # Initialize lists to store intersections and shared treated counts for treated unit `i`
#     intersected_controls <- list()
#     treated_with_shared_controls <- 0
#
#     for (j in 1:N1) {
#       if (i != j) {
#         # Find the shared controls between treated unit `i` and treated unit `j`
#         shared <- intersect(matched_pairs[[i]]$control_indices, matched_pairs[[j]]$control_indices)
#         intersected_controls[[j]] <- shared
#
#         # Count if there are any shared controls with treated unit `j`
#         if (length(shared) > 0) {
#           treated_with_shared_controls <- treated_with_shared_controls + 1
#         }
#       }
#     }
#
#     # Flatten the list of intersected controls and count unique values
#     unique_shared_controls <- unique(unlist(intersected_controls))
#     shared_controls[i] <- length(unique_shared_controls) # Store the count
#
#     # Store the number of treated units sharing controls with treated unit `i`
#     shared_treated[i] <- treated_with_shared_controls
#   }
#
#   list(
#     avg_shared_controls = mean(shared_controls),
#     p75_shared_controls = quantile(shared_controls, probs = 0.75),
#     avg_shared_treated = mean(shared_treated),
#     p75_shared_treated = quantile(shared_treated, probs = 0.75)
#   )
# }

generate_file_name <- function(
    N1 = 5,
    M = 10,
    i = 1,
    panel = "A"){
  file_name <- paste0("1d_DGP-", panel, "_M-", M, "_N1-", N1, "_i-", i, ".rds")
  return(file_name)
}

generate_one_dgp_and_matched_table <- function(N1 = 5,
                                        N0 = 1000,
                                        M = 10,
                                        i = 1,
                                        panel = "A",
                                        seed = NA,
                                        verbose = 0) {
  if (is.na(seed)){
    set.seed(123 + 12 * i)
  }

    if (verbose >= 1) {
      cat(sprintf("Generating DGP for replicate %d (Panel %s, M = %d, N1 = %d)...\n", i, panel, M, N1))
    }
    dgp_data <- generate_dgp(N1, N0, panel)
    file_name <- generate_file_name(N1 = N1,M = M,i = i,panel = panel)
    # file_name <- paste0("1d_DGP-", panel, "_M-", M, "_N1-", N1, "_i-", i, ".rds")
    if (verbose >= 2) {
      cat(sprintf("Saving matched table to %s\n", file_name))
    }
    generate_full_matched_table(dat = dgp_data,
                                M = M,
                                file_name = file_name)
}

read_one_matched_table <- function(N1 = 5,
                                   M = 10,
                                   i = 1,
                                   panel = "A"){
  file_name <- generate_file_name(N1 = N1,M = M,i = i,panel = panel)
  full_matched_table <- readRDS(file = here("scripts/new-inference/data/",file_name))
  return(full_matched_table)
}

generate_all_dgp_and_matched_table <- function(
    N0 = 1000,
    N1_values = c(5,10),
    M_values = c(1,4,10),
    panels = c("A","B"),
    num_replicates = 1000,
    verbose = 0){
  for (panel in panels) {
    for (N1 in N1_values) {
      for (M in M_values) {
        for (i in 1:num_replicates) {
          file_name <- generate_file_name(N1 = N1, M = M, i = i, panel = panel)
          if (file.exists(file = here("scripts/new-inference/data/",file_name))) {
            if (verbose == 2) {
              cat(sprintf("File '%s' exists, skipping.\n", file_name))
            }
          } else {
            generate_one_dgp_and_matched_table(
              N1 = N1,
              N0 = N0,
              M = M,
              i = i,
              panel = panel,
              verbose = verbose
            )
          }

        }
      }
    }
  }
}


compute_rejection_rate <-
  function(N1, N0, M, tau_0, alpha,
           num_replicates, max_permutations,
           panel,
           overlap_adj = T,
           sd_method = "usual") {
  rejection_rates <- numeric(num_replicates)
  overlap_stats <- list()

  for (i in 1:num_replicates) {
   full_matched_table <- read_one_matched_table(N1 = N1,
                                       M = M,
                                       i = i,
                                       panel = panel)
   rejection_rates[i] <-
     ferman_sign_change_test(full_matched_table = full_matched_table,
                             tau_0 = tau_0,
                             alpha = alpha,
                             max_permutations = max_permutations,
                             overlap_adj = overlap_adj,
                             sd_method = sd_method)

   matched_matrix <-
     get_matched_control_ids(full_matched_table)
   shared_neighbors <-
     compute_shared_neighbors(matched_matrix)
   overlap_stats[[i]] <-
     compute_overlap_statistics(shared_neighbors)
   # rejection_rates[i] <- ferman_sign_change_test(matched_pairs, dgp_data$treated, tau_0, alpha, max_permutations)
    # overlap_stats[[i]] <- compute_overlap_statistics(matched_pairs)
  }

  avg_rejection_rate <- mean(rejection_rates)

  # Aggregate overlap statistics
  avg_shared_controls <- mean(sapply(overlap_stats, `[[`, "avg_shared_controls"))
  p75_shared_controls <- mean(sapply(overlap_stats, `[[`, "p75_shared_controls"))
  avg_shared_treated <- mean(sapply(overlap_stats, `[[`, "avg_shared_treated"))
  p75_shared_treated <- mean(sapply(overlap_stats, `[[`, "p75_shared_treated"))

  list(
    rejection_rate = avg_rejection_rate,
    overlap_statistics = list(
      avg_shared_controls = avg_shared_controls,
      p75_shared_controls = p75_shared_controls,
      avg_shared_treated = avg_shared_treated,
      p75_shared_treated = p75_shared_treated
    )
  )
}


generate_full_table <- function(
    N0,
    N1_values,
    M_values,
    panels,
    tau_0,
    alpha,
    num_replicates,
    max_permutations,
    result_path = here("scripts/new-inference/outputs/results_table.rds")) {
  results <- list()

  for (panel in panels) {
    for (N1 in N1_values) {
      for (M in M_values) {
          cat("Running Panel", panel, "N1 =", N1, "M =", M, "\n")
          result <- compute_rejection_rate(N1, N0, M, tau_0, alpha, num_replicates, max_permutations, panel)
          results[[paste(panel, N1, M, sep = "_")]] <- result
      }
    }
  }

  saveRDS(results,
          result_path)
  results
}


#
# Function to load and parse results
load_results <- function(results_path, table_type) {
  # results_path <- here("scripts/new-inference/outputs/results_table-8-Dec.rds")
  results <- readRDS(results_path)

  # Extract relevant data based on the table type
  data <- lapply(results, function(x) {
    switch(
      table_type,
      "rejection_rate" = x$rejection_rate,
      "time_used" = x$time_used,
      "avg_shared_controls" = x$overlap_statistics$avg_shared_controls,
      "p75_shared_controls" = x$overlap_statistics$p75_shared_controls,
      "avg_shared_treated" = x$overlap_statistics$avg_shared_treated,
      "p75_shared_treated" = x$overlap_statistics$p75_shared_treated,
      stop("Invalid table type!")
    )
  })

  return(list(data = data, results = results))
}


# Function to reshape results into table format with panels and N1 rows
reshape_results_hierarchical <- function(data, results) {
  panels <- unique(sapply(names(results), function(name) strsplit(name, "_")[[1]][1]))
  n1_values <- unique(sapply(names(results), function(name) strsplit(name, "_")[[1]][2]))
  m_values <- unique(sapply(names(results), function(name) strsplit(name, "_")[[1]][3]))

  # Create a hierarchical structure
  table_data <- list()

  for (panel in panels) {
    panel_data <- matrix(NA, nrow = length(n1_values), ncol = length(m_values))
    colnames(panel_data) <- paste("M =", m_values)
    rownames(panel_data) <- paste("N1 =", n1_values)

    for (i in seq_along(n1_values)) {
      for (j in seq_along(m_values)) {
        result_key <- paste(panel, n1_values[i], m_values[j], sep = "_")
        if (result_key %in% names(data)) {
          panel_data[i, j] <- data[[result_key]]
        }
      }
    }
    table_data[[panel]] <- panel_data
  }

  return(table_data)
}

# Function to generate LaTeX code with super rows for panels
generate_latex_code_hierarchical <- function(table_data, caption, label) {
  latex_code <- paste0("\\begin{table}[H]\n",
                       "\\centering\n",
                       "\\caption{", caption, "}\n",
                       "\\label{", label, "}\n",
                       "\\begin{tabular}{l", paste(rep("c", ncol(table_data[[1]])), collapse = ""), "}\n",
                       "\\toprule\n")

  # Header row for M values
  latex_code <- paste0(latex_code, " & ", paste(colnames(table_data[[1]]), collapse = " & "), " \\\\\n", "\\midrule\n")

  # Add data for each panel
  for (panel_name in names(table_data)) {
    panel_matrix <- table_data[[panel_name]]
    latex_code <- paste0(latex_code, "\\multicolumn{", ncol(panel_matrix) + 1, "}{l}{\\textit{Panel ", panel_name, "}} \\\\\n")
    for (i in 1:nrow(panel_matrix)) {
      latex_code <- paste0(latex_code, rownames(panel_matrix)[i], " & ",
                           paste(round(panel_matrix[i, ], 3), collapse = " & "), " \\\\\n")
    }
    latex_code <- paste0(latex_code, "\\midrule\n")
  }

  latex_code <- paste0(latex_code, "\\bottomrule\n\\end{tabular}\n\\end{table}")
  return(latex_code)
}

# Updated Master function for LaTeX table
create_latex_table_hierarchical <- function(results_path, table_type, caption, label) {
  parsed_results <- load_results(results_path, table_type)
  table_data <- reshape_results_hierarchical(parsed_results$data, parsed_results$results)
  latex_code <- generate_latex_code_hierarchical(table_data, caption, label)
  return(latex_code)
}


# # Function to implement Ferman's Sign-changes algorithm
# ferman_sign_change_test_non_adjusted <- function(dgp_data, M, tau_0, alpha = 0.05, max_permutations = 1000) {
#   treated <- dgp_data$treated
#   control <- dgp_data$control
#
#   # Step 1: Match M nearest neighbors
#   treated <- treated %>%
#     rowwise() %>%
#     mutate(
#       matched_controls = list(
#         control %>%
#           arrange(abs(X - X[ID])) %>%
#           slice_head(n = M) %>%
#           pull(Y)
#       )
#     )
#
#   # Step 2: Compute individual treatment effects
#   treated <- treated %>%
#     mutate(
#       tau_hat = Y - mean(unlist(matched_controls)) - tau_0
#     )
#
#   tau_hat_S <- mean(treated$tau_hat)
#
#   # Step 3: Compute test statistic
#   T_obs <- abs(tau_hat_S) / sd(treated$tau_hat)
#
#   # Step 4: Null distribution
#   S <- treated$tau_hat
#   G <- replicate(max_permutations, sample(c(-1, 1), length(S), replace = TRUE))
#   T_null <- apply(G, 2, function(g) {
#     gS <- g * S
#     mean_gS <- mean(gS)
#     abs(mean_gS) / sqrt(sum((gS - mean_gS)^2) / (length(S) - 1))
#   })
#
#   # Step 5: Decision rule
#   critical_value <- quantile(T_null, probs = 1 - alpha)
#   reject <- as.numeric(T_obs > critical_value)
#
#   reject
# }
