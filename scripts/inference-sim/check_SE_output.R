n_t_values <- c(5, 10, 25, 50)
beta_c_values <- c(0, 2, 4)
se_estimates <-
  readRDS(here::here("data/outputs/inference-sim",
                     "se_estimates_check_SE.rds"))

adjusted_se_estimates <- se_estimates
for (i in seq_along(n_t_values)) {
  for (j in seq_along(beta_c_values)) {
    n_t <- n_t_values[i]

    # Adjust SE estimates for Boot SE (Sign, Uniform Weight) and Boot SE (Sign, SCM Weight)
    for (k in c(3, 4)) {
      adjusted_se_estimates[i, j, , k] <- se_estimates[i, j, , k] / sqrt(n_t)
    }
  }
}


for (i in seq_along(n_t_values)) {
  for (j in seq_along(beta_c_values)) {
    for (k in 1:6) {
      if (k == 1){
        bias_matrix[i, j, k] <- mean(se_estimates[i, j, , 1])
      }else{
        bias_matrix[i, j, k] <- mean(adjusted_se_estimates[i, j, , k] - se_estimates[i, j, , 1])
      }

    }
  }
}

saveRDS(bias_matrix,
        here::here("data/outputs/inference-sim",
                   "bias_matrix_check_SE.rds"))

bias_matrix <-
  readRDS(here::here("data/outputs/inference-sim",
                     "bias_matrix_check_SE.rds"))

estimator_names <- c(
  "True CMSE",
  "AE SE",
  "Boot SE (Sign, Uniform Weight)",
  "Boot SE (Sign, SCM Weight)",
  "Boot SE (Wild, Uniform Weight)",
  "Boot SE (Wild, SCM Weight)"
)

# Change the column names to "High Overlap," "Mid Overlap," and "Low Overlap"
colnames(bias_matrix) <- c("High Overlap", "Mid Overlap", "Low Overlap")

# Function to generate LaTeX code for all bias matrices combined into one table
generate_combined_latex_table <- function(bias_matrix, estimator_names, n_t_values, beta_c_values) {
  latex_content <- ""

  latex_content <- paste0(latex_content, "\\begin{table}[h]\n")
  latex_content <- paste0(latex_content, "\\centering\n")
  latex_content <- paste0(latex_content, "\\begin{tabular}{|c|", paste(rep("c", length(beta_c_values)), collapse = "|"), "|}\n")
  latex_content <- paste0(latex_content, "\\hline\n")

  for (k in 1:6) {  # Skip the first estimator ("True CMSE") by starting at 2
    estimator_name <- estimator_names[k]  # Adjust index since we start at 2
    bias_matrix_est <- bias_matrix[, , k]

    latex_content <- paste0(latex_content, "\\multicolumn{", length(beta_c_values) + 1, "}{|c|}{", estimator_name, "} \\\\\n")
    latex_content <- paste0(latex_content, "\\hline\n")
    latex_content <- paste0(latex_content, "& ", paste(colnames(bias_matrix_est), collapse = " & "), "\\\\\n")
    latex_content <- paste0(latex_content, "\\hline\n")

    for (i in 1:length(n_t_values)) {
      row_name <- gsub("_", "\\\\_", rownames(bias_matrix_est)[i])  # Replace _ with \_
      latex_content <- paste0(latex_content, row_name, " & ", paste(signif(bias_matrix_est[i, ], 2), collapse = " & "), "\\\\\n")
    }
    latex_content <- paste0(latex_content, "\\hline\n")
  }

  latex_content <- paste0(latex_content, "\\end{tabular}\n")
  latex_content <- paste0(latex_content, "\\caption{Bias Matrix for Different Estimators}\n")
  latex_content <- paste0(latex_content, "\\end{table}\n\n")

  return(latex_content)
}

# Generate LaTeX content
latex_table <- generate_combined_latex_table(bias_matrix, estimator_names, n_t_values, colnames(bias_matrix))

# Write LaTeX table to a .txt file
output_file <- here::here("scripts/inference-sim", "bias_matrix_check_SE.txt")
writeLines(latex_table, output_file)

