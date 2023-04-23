
# note: standard, bayesian, and wild bootstraps all capture true sd
#  - but wild bootstrap requires vec-mean(vec) 
#    BEFORE multiplying by the weight
#  - others work with vec-mean(vec), but it's not required

N <- 100
B <- 1000
mu <- 1

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

# bayesian bootstrap

map_dbl(
  1:B,
  function(x) {
    # (vec * as.numeric(gtools::rdirichlet(1, alpha=rep(1,N)))) %>% sum()
    (vec * as.numeric(rmultinom(1, size=N, prob=rep(1/N,N))) / N) %>% sum()
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

