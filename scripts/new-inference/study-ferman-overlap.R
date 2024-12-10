library(here)
source(here("scripts/new-inference/utils-replicate-ferman.R"))

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
