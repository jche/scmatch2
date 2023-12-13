
# note: standard, bayesian, and wild bootstraps all capture true sd
#  - but wild bootstrap requires vec-mean(vec) 
#    BEFORE multiplying by the weight
#  - others work with vec-mean(vec), but it's not required
library(tidyverse)
N <- 100
B <- 1000
mu <- 1



I <- 100 # do 100 times
covered <- CI_lower <- CI_upper <- numeric(I)
T_star <- numeric(B)
set.seed(123)
for (i in 1:I){
  print(i)
  # Generate data of size 100 according to the DGP N(mu,1)
  X <- rnorm(N, mean=mu); X_bar = mean(X)
  # Get the residuals 
  X_resid = X - X_bar
  # Perform Bayesian bootstrap
  for (b in 1:B){
    W = gtools::rdirichlet(1, alpha=rep(1,N))
    T_star[b] = sum(X_resid * W)
  }
  
  CI_lower[i] = X_bar - quantile(T_star, 0.975)
  CI_upper[i] = X_bar - quantile(T_star, 0.025)
  covered[i] = (CI_lower[i] < mu) & (mu < CI_upper[i])
}
mean(covered)
idx_not_covered<-which(covered!=T)
CI_lower[idx_not_covered]
CI_upper[idx_not_covered]


# goal: estimate mean
vec <- rnorm(N, mean=mu)

# true sd:
map_dbl(
  1:B,
  function(x) {
    rnorm(N) %>% mean(mean=mu)
  }
) %>% 
  sd()

# naive bootstrap
map_dbl(
  1:B,
  function(x) {
    sample(vec, N, replace=T) %>% mean()
  }) %>% 
  sd()

# non-par bootstrap, which is equivalent to naive bootstrap
map_dbl(
  1:B,
  function(x) {
    (vec * as.numeric(rmultinom(1, size=N, prob=rep(1/N,N))) / N) %>% sum()
  }) %>% 
  sd()

# bayesian bootstrap (Rubin, 1981)
map_dbl(
  1:B,
  function(x) {
    (vec * as.numeric(gtools::rdirichlet(1, alpha=rep(1,N)))) %>% sum()
  }) %>% 
  sd()

# wild bootstrap
#  - note: needs to subtract mean to be accurate!
map_dbl(
  1:B,
  function(x) {
    ((vec-mean(vec)) * sample(
      c( -(sqrt(5)-1)/2, (sqrt(5)+1)/2 ),
      prob = c( (sqrt(5)+1)/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5)) ),
      replace = T, size = N)) %>% mean()
  }) %>% 
  sd()

