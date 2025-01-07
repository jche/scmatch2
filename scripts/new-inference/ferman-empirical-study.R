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
## Toy example of sign test.
#######
hist(X)

T_null <- replicate(1000, {
  signs <- sample(c(-1, 1), N1+N0, replace=T)
  var(signs * X)
})
hist(T_null)
var(X)


######
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
