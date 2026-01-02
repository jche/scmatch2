library(here)
source(here("scripts/new-inference/utils-replicate-ferman.R"))

######
### Verification that complex epsilon_1 is not symmetrical
#####
epsilon_1 <- (rchisq(1000, df = 1) - 1) / sqrt(2)
hist(epsilon_1)
abs( mean(epsilon_1) - 0 ) <= 0.01

######
### Verification that complex mu_1 is not symmetrical
#####
N1 = 10; N0 = 1000
X <- c(rnorm(N1, mean = 0, sd = 1),
        rnorm(N0, mean = 0, sd = 1))
mu_1 <- function(x) (qchisq(pnorm(x), df = 8) - 8) / sqrt(16)
hist(mu_1(X))
mean(mu_1(X))

########
## Showing that variance in the signed test is random
#######
N1 = 10; N0 = 1000
X <- c(rnorm(N1, mean = 0, sd = 1),
       rnorm(N0, mean = 0, sd = 1))
hist(X)

var_dist_X_fixed <- replicate(1000, {
  signs <- sample(c(-1, 1), N1+N0, replace=T)
  var(signs * X)
})
hist(var_dist_X_fixed)
var(var_dist_X_fixed) # 1.456992e-06

var_dist_X_random <- replicate(1000, {
  X <- c(rnorm(N1, mean = 0, sd = 1),
         rnorm(N0, mean = 0, sd = 1))
  signs <- sample(c(-1, 1), N1+N0, replace=T)
  var(signs * X)
})
hist(var_dist_X_random)
var(var_dist_X_random) # 0.001960424
var(var_dist_X_random) / var(var_dist_X_fixed) # 1345 times


######
## Comparison of methods
#######
######
## Check the method difference
######
## Goal: save two matrices of (num_replicates x max_permutations) of variance and
##  from sd_method == "usual" and
## With that we can later compute the critical value
## The critical value is irrelavant of how the observed
## There are two ways of doing it: a) make it a part of; b) or write directly
## Decision: choose b) because a) requires too much testing
N1 = 10
M = 4
panel = "D"
num_replicates = 300
max_permutations = 200
library(here)
source(here("scripts/new-inference/utils-replicate-ferman.R"))
load_all()


col_names <- c("iteration", "sd_method", paste0("null_dist_", 1:max_permutations))
null_dist_df <- data.frame(matrix(NA, nrow = 2 * num_replicates, ncol = length(col_names)))
colnames(null_dist_df) <- col_names


save_interval <- 10
output_data_path <- here("scripts/new-inference/outputs/null_dist_df.rds")
output_time_path <- here("scripts/new-inference/outputs/null_dist_df_computation_time.rds")

# Initialize or load existing data
if (file.exists(output_data_path)) {
  null_dist_df <- readRDS(output_data_path)
  row_index <- max(which(!is.na(null_dist_df$iteration))) + 1  # Find the next empty row
} else {
  col_names <- c("iteration", "sd_method", paste0("null_dist_", 1:max_permutations))
  null_dist_df <- data.frame(matrix(NA, nrow = 2 * num_replicates, ncol = length(col_names)))
  colnames(null_dist_df) <- col_names
  row_index <- 1  # Row index tracker
}


tstart <- proc.time()
for (i in ceiling(row_index / 2):num_replicates) {
  # i <- 1
  set.seed(123 + 12 * i)
  full_matched_table <- read_one_matched_table(N1 = N1,
                                               M = M,
                                               i = i,
                                               panel = panel)

  T_null_original <- ferman_sign_change_null_dist(
    full_matched_table = full_matched_table,
    tau_0 = 0,
    max_permutations = max_permutations,
    overlap_adj = T,
    sd_method = "usual")

  T_null_pooled <- ferman_sign_change_null_dist(
    full_matched_table = full_matched_table,
    tau_0 = 0,
    max_permutations = max_permutations,
    overlap_adj = T,
    sd_method = "pooled")

  # Assign results for "usual"
  null_dist_df[row_index, "iteration"] <- i
  null_dist_df[row_index, "sd_method"] <- "usual"
  null_dist_df[row_index, paste0("null_dist_", 1:max_permutations)] <- T_null_original
  row_index <- row_index + 1

  # Assign results for "pooled"
  null_dist_df[row_index, "iteration"] <- i
  null_dist_df[row_index, "sd_method"] <- "pooled"
  null_dist_df[row_index, paste0("null_dist_", 1:max_permutations)] <- T_null_pooled
  row_index <- row_index + 1

  # Verbose printing
  if (i %% save_interval == 0 || i == num_replicates) {
    cat(sprintf("Progress: Completed replicate %d of %d\n", i, num_replicates))
    saveRDS(null_dist_df, output_data_path)  # Save intermediate results
  }
}

tend <- proc.time()

# Save final outputs
saveRDS(null_dist_df, output_data_path)
saveRDS(tend - tstart, output_time_path)

cat("Final results saved.\n")

#######
# The original comparison code
#######
result_ferman <-
  compute_rejection_rate(
    N1, N0, M, tau_0, alpha,
    num_replicates=300,
    max_permutations = 200,
    panel,
    overlap_adj = T)
source(here("scripts/new-inference/utils-replicate-ferman.R"))
result_pooled <-
  compute_rejection_rate(
    N1, N0, M, tau_0, alpha,
    num_replicates=300,
    max_permutations = 200,
    panel,
    overlap_adj = T,
    sd_method = "pooled")
# This runs 1 hour
saveRDS(result_ferman, here("scripts/new-inference/outputs/result_ferman.rds"))
saveRDS(result_pooled, here("scripts/new-inference/outputs/result_pooled.rds"))

result_pooled <- readRDS(here::here("scripts/new-inference/outputs/result_pooled.rds"))
result_ferman <- readRDS(here::here("scripts/new-inference/outputs/result_ferman.rds"))





#####
## Study of overlap
####
######
####
## Check the overlap impact
N1 = 10
M = 4
panel = "D"
cat("Running Panel", panel, "N1 =", N1, "M =", M, "\n")
result_has_overlap_adj <-
  compute_rejection_rate(
    N1, N0, M, tau_0, alpha,
    num_replicates=300,
    max_permutations = 200,
    panel)
source(here("scripts/new-inference/utils-replicate-ferman.R"))
result_no_overlap_adj <-
  compute_rejection_rate(
    N1, N0, M, tau_0, alpha,
    num_replicates=300,
    max_permutations=200,
    panel,
    overlap_adj = F)
## Seems no impact -- i think this needs to change to N1 = 50, M = 10
##  where overlap is more significant



## Under null, is the null distribution of tau(S)
## and \tau(gS) the same or different?
#######
full_matched_table <-
  read_one_matched_table(N1 = 50,
                         M = 10,
                         i = 1,
                         panel = "E")
null_dist_adj <-
  ferman_sign_change_null_dist(full_matched_table,
                               tau_0 = 0,
                               max_permutations = 1000,
                               overlap_adj = T)
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

#####
### Study of overlap
full_matched_table <-
  read_one_matched_table(N1 = 50,
                         M = 10,
                         i = 1,
                         panel = "E")
## input the full_matched_table, output the distribution of the
source(here("scripts/new-inference/utils-replicate-ferman.R"))
null_dist_adj <-
  ferman_sign_change_null_dist(full_matched_table,
           tau_0 = 0,
           max_permutations = 1000,
           overlap_adj = T)
null_dist_no_adj <-
  ferman_sign_change_null_dist(full_matched_table,
                               tau_0 = 0,
                               max_permutations = 1000,
                               overlap_adj = F)

par(mfrow = c(1,2))
hist(null_dist_adj)
hist(null_dist_no_adj)
par(mfrow = c(1,1))

alpha = 0.1
quantile(null_dist_adj, probs = 1 - alpha)
quantile(null_dist_no_adj, probs = 1 - alpha)
## this says there is not a big difference with and
## without adjustment

# Generate valid sign changes
sum_signs_adj <-function(overlap_adj){
  replicate(max_permutations, {
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
    sum(signs)
  })
}

signs_adj <- sum_signs_adj(T)
signs_no_adj <- sum_signs_adj(F)
summary(signs_adj)
summary(signs_no_adj)

## Let's see overlap proportions
# Basically the overlap issue is little
matched_matrix <- get_matched_control_ids(full_matched_table)
shared_neighbors <- compute_shared_neighbors(matched_matrix)
shared_neighbors_binary <- shared_neighbors > 0
compute_overlap_statistics(shared_neighbors)
