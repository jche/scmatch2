library(MASS) # For mvrnorm
library(dplyr) # For data manipulation
library(here)

# Function to generate DGP based on the panel
generate_dgp <- function(N1, N0, panel = "A") {
  # Generate covariates
  X_treated <- rnorm(N1, mean = 0, sd = 1)
  X_control <- rnorm(N0, mean = 0, sd = 1)

  # Generate outcomes based on the panel
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
  }

  Y_treated <- sapply(X_treated, function(x) mu_1(x) + epsilon_1())
  Y_control <- rnorm(N0, mean = 0, sd = 1)

  list(
    treated = data.frame(ID = 1:N1, X = X_treated, Y = Y_treated),
    control = data.frame(ID = 1:N0, X = X_control, Y = Y_control)
  )
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

match_controls <- function(treated, control, M) {
  matched_pairs <- lapply(1:nrow(treated), function(i) {
    distances <- abs(control$X - treated$X[i])
    matched_indices <- order(distances)[1:M]
    list(
      controls = control$Y[matched_indices],
      control_indices = matched_indices
    )
  })
  return(matched_pairs)
}


ferman_sign_change_test <- function(matched_pairs, treated, tau_0, alpha = 0.05, max_permutations = 1000) {
  N1 <- nrow(treated)

  # Calculate individual treatment effects
  tau_i <- sapply(1:N1, function(i) {
    treated$Y[i] - mean(matched_pairs[[i]]$controls) - tau_0
  })

  tau_bar <- mean(tau_i)

  # Test statistic
  T_obs <- abs(tau_bar) / sqrt(sum((tau_i - tau_bar)^2) / (N1 - 1))

  # Generate sign changes respecting shared neighbors
  shared_neighbors <- matrix(FALSE, N1, N1)
  for (i in 1:N1) {
    for (j in 1:N1) {
      if (i < j) {
        shared <- intersect(matched_pairs[[i]]$control_indices, matched_pairs[[j]]$control_indices)
        shared_neighbors[i, j] <- shared_neighbors[j, i] <- length(shared) > 0
      }
    }
  }

  # Generate valid sign changes
  T_null <- replicate(max_permutations, {
    signs <- sample(c(-1, 1), 1)
    for (i in 2:N1) {
      connected <- which(shared_neighbors[i, 1:(i - 1)])
      if (length(connected) > 0) {
        signs[i] <- signs[connected[1]]
      } else {
        signs[i] <- sample(c(-1, 1), 1)
      }
    }
    gS <- signs * tau_i
    gS_bar <- mean(gS)
    abs(gS_bar) / sqrt(sum((gS - gS_bar)^2) / (N1 - 1))
  })

  critical_value <- quantile(T_null, probs = 1 - alpha)
  reject <- as.numeric(T_obs > critical_value)
  reject
}


# Compute overlap-related statistics
compute_overlap_statistics <- function(matched_pairs) {
  N1 <- length(matched_pairs)

  shared_controls <- numeric(N1) # Initialize a vector to store shared control counts
  shared_treated <- numeric(N1)  # Initialize a vector to store shared treated counts

  for (i in 1:N1) {
    # Initialize lists to store intersections and shared treated counts for treated unit `i`
    intersected_controls <- list()
    treated_with_shared_controls <- 0

    for (j in 1:N1) {
      if (i != j) {
        # Find the shared controls between treated unit `i` and treated unit `j`
        shared <- intersect(matched_pairs[[i]]$control_indices, matched_pairs[[j]]$control_indices)
        intersected_controls[[j]] <- shared

        # Count if there are any shared controls with treated unit `j`
        if (length(shared) > 0) {
          treated_with_shared_controls <- treated_with_shared_controls + 1
        }
      }
    }

    # Flatten the list of intersected controls and count unique values
    unique_shared_controls <- unique(unlist(intersected_controls))
    shared_controls[i] <- length(unique_shared_controls) # Store the count

    # Store the number of treated units sharing controls with treated unit `i`
    shared_treated[i] <- treated_with_shared_controls
  }

  list(
    avg_shared_controls = mean(shared_controls),
    p75_shared_controls = quantile(shared_controls, probs = 0.75),
    avg_shared_treated = mean(shared_treated),
    p75_shared_treated = quantile(shared_treated, probs = 0.75)
  )
}

compute_rejection_rate <- function(N1, N0, M, tau_0, alpha, num_replicates, max_permutations, panel) {
  rejection_rates <- numeric(num_replicates)
  overlap_stats <- list()

  for (i in 1:num_replicates) {
    dgp_data <- generate_dgp(N1, N0, panel)
    matched_pairs <- match_controls(dgp_data$treated, dgp_data$control, M)
    rejection_rates[i] <- ferman_sign_change_test(matched_pairs, dgp_data$treated, tau_0, alpha, max_permutations)
    overlap_stats[[i]] <- compute_overlap_statistics(matched_pairs)
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


# Example of full table generation
generate_full_table <- function(
    N0, N1_values, M_values, panels, tau_0, alpha, num_replicates, max_permutations,
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

# # Function to reshape results into a table format
# reshape_results <- function(data, results) {
#   panels <- unique(sapply(names(results), function(name) strsplit(name, "_")[[1]][1]))
#   n1_values <- unique(sapply(names(results), function(name) strsplit(name, "_")[[1]][2]))
#   m_values <- unique(sapply(names(results), function(name) strsplit(name, "_")[[1]][3]))
#
#   table_data <- matrix(NA, nrow = length(n1_values), ncol = length(m_values) * length(panels))
#   col_names <- vector()
#
#   col_idx <- 1
#   for (panel in panels) {
#     for (m in m_values) {
#       col_names <- c(col_names, paste("Panel", panel, "M =", m))
#       for (i in seq_along(n1_values)) {
#         result_key <- paste(panel, n1_values[i], m, sep = "_")
#         table_data[i, col_idx] <- ifelse(result_key %in% names(data), data[[result_key]], NA)
#       }
#       col_idx <- col_idx + 1
#     }
#   }
#
#   row_names <- paste("N1 =", n1_values)
#   table_data <- as.data.frame(table_data)
#   names(table_data) <- col_names
#   rownames(table_data) <- row_names
#
#   return(table_data)
# }
#
# # Function to generate LaTeX table
# generate_latex_code <- function(table_data, caption, label) {
#   latex_code <- paste0("\\begin{table}[H]\n",
#                        "\\centering\n",
#                        "\\caption{", caption, "}\n",
#                        "\\label{", label, "}\n",
#                        "\\begin{tabular}{l", paste(rep("c", ncol(table_data)), collapse = ""), "}\n",
#                        "\\hline\n",
#                        " & ", paste(names(table_data), collapse = " & "), " \\\\\n",
#                        "\\hline\n")
#
#   for (i in seq_len(nrow(table_data))) {
#     latex_code <- paste0(latex_code, rownames(table_data)[i], " & ",
#                          paste(round(table_data[i, ], 3), collapse = " & "), " \\\\\n")
#   }
#
#   latex_code <- paste0(latex_code, "\\hline\n\\end{tabular}\n\\end{table}")
#   return(latex_code)
# }
#
# # Master function to create LaTeX table
# create_latex_table <- function(results_path, table_type, caption, label) {
#   parsed_results <- load_results(results_path, table_type)
#   table_data <- reshape_results(parsed_results$data, parsed_results$results)
#   latex_code <- generate_latex_code(table_data, caption, label)
#   return(latex_code)
# }
#

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
    latex_code <- paste0(latex_code, "\\multicolumn{", ncol(panel_matrix) + 1, "}{l}{\\textit{", panel_name, "}} \\\\\n")
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
