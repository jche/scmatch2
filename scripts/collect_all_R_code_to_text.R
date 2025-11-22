# Set the folder path (change to your desired directory if needed)
folder_path <- here::here("scripts/ferman-analysis")

# Get list of all .R files in the folder (non-recursive)
r_files <- list.files(path = folder_path,
                      pattern = "\\.R$",
                      full.names = TRUE)

# # Define the subset of R files you want in the order you want them
# r_files <- file.path(folder_path, c(
#   "01_fit_models.R",
#   "02_compute_hazards.R",
#   "03_dr_estimator.R",
#   "04_compute_rmst_eta.R",
#   "05_dr_pipeline.R",
#   "datagen.R",
#   "utils.R"
# ))



# Initialize output file
output_file <- file.path(folder_path, "all_scripts_combined.txt")
file.create(output_file)

# Loop through each R file and append its contents to the output file
for (r_file in r_files) {
  cat("\n### START OF FILE:", basename(r_file), "###\n", file = output_file, append = TRUE)
  file_content <- readLines(r_file, warn = FALSE)
  cat(file_content, sep = "\n", file = output_file, append = TRUE)
  cat("\n### END OF FILE:", basename(r_file), "###\n\n", file = output_file, append = TRUE)
}

cat("✅ All R scripts have been collected into:", output_file, "\n")
