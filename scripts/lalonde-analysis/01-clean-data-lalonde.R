# scripts/lalonde-analysis/01-clean-data-lalonde.R
# Prep raw data for analysis by renaming variables to match example-lalonde.R logic

library(tidyverse)
library(here)

# Load data based on the example script preference
# Note: Adjust path if your project root differs, assuming standard structure
load(file = here::here("data/inputs/lalonde_w_cps_cos.RData"))

# In the loaded RData, the dataframe is usually named 'lalonde_df'
# We apply the renaming logic from example-lalonde.R
lalonde_clean <- lalonde_df %>%
  rename(
    X1 = black,
    X2 = hispanic,
    X3 = married,
    X4 = nodegree,
    X5 = age,
    X6 = education,
    X7 = re74,
    X8 = re75
  ) %>%
  mutate(
    Y = Y, # Assuming Y is the outcome
    Z = Z     # Treatment is already Z
  ) %>%
  select(Z, Y, starts_with("X"))

# Create directory if it doesn't exist
if(!dir.exists(here::here("scripts/lalonde-analysis/data"))) {
  dir.create(here::here("scripts/lalonde-analysis/data"), recursive = TRUE)
}

# Save for the core analysis step
saveRDS(lalonde_clean,
        file = here::here("scripts/lalonde-analysis/data/lalonde_for_analysis.rds"))


