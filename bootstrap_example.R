
# simple example to unpack what's going on 
# with the bootstrap for the att
# - with no selection, random treatment effect

# Initial conclusions:
#  - regular boot & multinomial boot differ because 
#    regular boot is using bootstrapped sum(Z), 
#    while multinomial boot implicitly conditions 
#    on n1=sum(Z) in the original data.


require(tidyverse)

# generate data
numsamp <- 1000
tau <- 1
df <- tibble(
  Y0 = rnorm(numsamp),
  tauj = rnorm(numsamp, mean=tau),
  Y1 = Y0 + tauj,
  Z = sample(0:1, size=numsamp, replace=T, prob=c(0.9,0.1)),
  Yobs = ifelse(Z, Y1, Y0)) %>% 
  # mutate(wt = ifelse(Z, 1/sum(Z), 1/(n()-sum(Z))))
  mutate(wt = ifelse(Z, 1, sum(Z)/(n()-sum(Z))))

tau_hat <- df %>% 
  summarize(att = sum((2*Z-1) * wt * Yobs) / sum(Z))
tau_hat

# ground truth ------------------------------------------------------------

# true sd: randomness only in Z
#  - weird w.r.t. the ATT in some sense...
B <- 500
taus_resampZ <- map_dbl(
  1:B,
  function(x) {
    df %>% 
      mutate(Z = sample(0:1, size=numsamp, replace=T, prob=c(0.9,0.1)),
             Yobs = ifelse(Z, Y1, Y0),
             wt = ifelse(Z, 1, sum(Z)/(n()-sum(Z)))) %>% 
      summarize(att = sum((2*Z-1) * wt * Yobs) / sum(Z)) %>% 
      pull(att)
  })
ggplot(tibble(x=taus_resampZ)) +
  geom_histogram(aes(x)) +
  geom_vline(xintercept=tau)
mean(taus_resampZ)
sd(taus_resampZ)

# true sd: randomness in Y
taus_resampY <- map_dbl(
  1:B,
  function(x) {
    df %>% 
      mutate(Y0 = rnorm(numsamp),
             tauj = rnorm(numsamp, mean=tau),
             Y1 = Y0 + tauj,
             Yobs = ifelse(Z, Y1, Y0),
             wt = ifelse(Z, 1, sum(Z)/(n()-sum(Z)))) %>% 
      summarize(att = sum((2*Z-1) * wt * Yobs) / sum(Z)) %>% 
      pull(att)
  })
ggplot(tibble(x=taus_resampY)) +
  geom_histogram(aes(x)) +
  geom_vline(xintercept=tau)
mean(taus_resampY)
sd(taus_resampY)

# note: given original dataset may not be similar to the "average" dataset!



# TODO: what is the "true" variation we're after?
#  - in some sense, we'd like to condition on the treated units...?
# regardless, bootstrap seems to be much closer to resampZ se than to resampY se...




# bootstrap results -------------------------------------------------------

# TODO CHECK: does it matter what the co units' weights sum to...?

boot <- function(d, B=100) {
  map_dbl(1:B,
          .progress = "Bootstrapping...",
          function(b) {
            d %>% 
              sample_n(n(), replace=T) %>%
              summarize(att = sum((2*Z-1) * wt * Yobs) / sum(Z)) %>%
              pull(att)
          })
}

boot2 <- function(d, B=100) {
  # n <- nrow(d)
  # n1 <- sum(d$Z)
  
  map_dbl(1:B,
          .progress = "Bootstrapping...",
          function(b) {
            att_hat <- d %>%
              summarize(att = sum((2*Z-1) * wt * Yobs) / sum(Z)) %>%
              pull(att)
            
            d %>% 
              sample_n(n(), replace=T) %>%
              summarize(att = sum((2*Z-1) * wt * Yobs - sum(Z)/n()*att_hat) / sum(Z)) %>%
              # summarize(att = sum((2*Z-1) * wt * Yobs - n1/n*att_hat) / n1) %>%
              pull(att)
          })
}

# use n1/n*att_hat to cancel bias for multinomial weights,
#  which sum to n instead of n1
#  - this actually helps the wild bootstrap, for whatever reason???
# NOTE: this is conservative, due to n1 weirdness
weighted_boot <- function(d, wt_fun, B=100) {
  n <- nrow(d)
  n1 <- sum(d$Z)
  
  att_hat <- d %>%
    summarize(att = sum((2*Z-1) * wt * Yobs) / n1) %>%
    pull(att)
  
  map_dbl(1:B,
          .progress = "Bootstrapping...",
          function(b) {
            d %>% 
              mutate(
                weights_boot = wt_fun(n),
                tau_i = (2*Z-1) * wt * Yobs) %>% 
              summarize(att = sum(weights_boot * 
                                    (tau_i - n1/n*att_hat)) / n1) %>%
              pull(att)
          })
}

# this is pretty close to what the bootstrap does!
weighted_boot2 <- function(d, wt_fun, B=100) {
  att_hat <- d %>%
    summarize(att = sum((2*Z-1) * wt * Yobs) / sum(Z)) %>%
    pull(att)
  
  map_dbl(1:B,
          .progress = "Bootstrapping...",
          function(b) {
            d %>% 
              mutate(
                weights_boot = wt_fun(n()),
                tau_i = (2*Z-1) * wt * Yobs) %>% 
              summarize(
                # n1 = sum(weights_boot*Z),
                # n = sum(weights_boot),
                att = sum(weights_boot * 
                            (tau_i - sum(weights_boot*Z)/n()*att_hat)) / 
                  sum(weights_boot*Z) ) %>%
              pull(att)
          })
}

# weighted_boot3 <- function(d, wt_fun, B=100) {
#   map_dbl(1:B,
#           .progress = "Bootstrapping...",
#           function(b) {
#             d %>% 
#               mutate(
#                 weights_boot = wt_fun(n()),
#                 tau_i = (2*Z-1) * wt * Yobs) %>% 
#               summarize(
#                 att = sum(weights_boot * (tau_i)) / sum(weights_boot*Z) ) %>%
#               pull(att)
#           })
# }

if (F) {
  # full bootstrap works great
  boot_samps <- boot(df, B=500)
  hist(boot_samps)
  sd(boot_samps)
  
  boot_samps2 <- boot2(df, B=500)
  hist(boot_samps2)
  sd(boot_samps2)

  # weighted bootstrap works...?
  
  wf <- function(n) as.numeric(rmultinom(1, size=n, prob=rep(1/n,n)))
  # wf <- function(n) as.numeric(gtools::rdirichlet(1, alpha=rep(1,n)))
  wf <- function(n) {
    sample(
      c( -(sqrt(5)-1)/2, (sqrt(5)+1)/2 ),
      prob = c( (sqrt(5)+1)/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5)) ),
      replace = T, size = n)
  }
  
  wboot_samps <- weighted_boot(df, wf, B=500)
  hist(wboot_samps)
  sd(wboot_samps)
  
  wboot2_samps <- weighted_boot2(df, wf, B=500)
  hist(wboot2_samps)
  sd(wboot2_samps)
  
  wboot3_samps <- weighted_boot3(df, wf, B=500)
  hist(wboot3_samps)
  sd(wboot3_samps)
  
  sd(boot_samps)
  sd(boot_samps2)
  sd(wboot_samps)
  sd(wboot2_samps)
  sd(wboot3_samps)
  sd(taus_resampY)
  sd(taus_resampZ)
}



