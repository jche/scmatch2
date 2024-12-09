library(FNN)
library(CSM)
library(dplyr)
library(here)
source(here("scripts/new-inference/utils-replicate-ferman.R"))

set.seed(4044440)

dat <- generate_dgp(N1 = 5, N0 = 1000, panel = "A")

### Using CSM::get_cal_matches
k = 4
# devtools::load_all()
mtch <- get_cal_matches(
  dat,
  covs = "X",
  treatment = "Z",
  scaling = 1,
  metric = "euclidean",
  rad_method = "knn",
  k = k
)

full_table <- full_unit_table(mtch,
                              nonzero_weight_only = F )
saveRDS(mtch, file = here("scripts/new-inference/data/one-full-table.rds"))




### Using FNN to do the matching
# Separate treated and control groups
treated <- dat %>% filter(Z == TRUE)
control <- dat %>% filter(Z == FALSE)

# Covariates for matching
treated_cov <- treated %>% dplyr::select(X)
control_cov <- control %>% dplyr::select(X)

# Perform KNN matching
k <- 4
knn_result <- get.knnx(data = control_cov,
                       query = treated_cov,
                       k = k)

# Extract indices of matches and compute matched control outcomes
matched_indices <- knn_result$nn.index # so this is not id, but the row number of the control data
matched_control_outcomes <- control$Y0[matched_indices]

# Create a result data frame
matching_result <- treated %>%
  mutate(
    matched_Y0 = rowMeans(matched_control_outcomes), # Average if k > 1
    ATT = Y - matched_Y0 # Compute ATT for each treated unit
  )

# Output results
print(matching_result)
