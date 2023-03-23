
# functions for simulating data


# toy example -------------------------------------------------------------

# generate data mostly in grid (0,1) x (0,1)
gen_df_adv <- function(nc, nt, 
                       f0_sd = 0.1,   # homoskedastic noise
                       f0_fun = function(X1, X2) {1},
                       tx_effect_fun = function(X1, X2) {1}
                       # f0_fun = function(X1, X2) { abs(X1-X2) },
                       # tx_effect = function(X1, X2) {(X1-0.5)^2+(X2-0.5)^2},
                       ) {
  
  SD <- 0.1
  
  # tx units clustered at (0.25,0.25) and (0.75,0.75)
  dat_txblobs <- tibble(
    X1 = c(rnorm(nt/2, mean=0.25, sd=SD), 
           rnorm(nt/2, mean=0.75, sd=SD)),
    X2 = c(rnorm(nt/2, mean=0.25, sd=SD), 
           rnorm(nt/2, mean=0.75, sd=SD)),
    Z  = T
  )
  
  # co units clustered at (0.75,0.25) and (0.25,0.75)
  dat_coblobs <- tibble(
    X1 = c(rnorm((nc-nt)/2, mean=0.25, sd=SD), 
           rnorm((nc-nt)/2, mean=0.75, sd=SD)),
    X2 = c(rnorm((nc-nt)/2, mean=0.75, sd=SD), 
           rnorm((nc-nt)/2, mean=0.25, sd=SD)),
    Z  = F
  )
  
  # some co units uniformly scattered on (0,1) box
  dat_conear <- tibble(
    X1 = runif(nt),
    X2 = runif(nt),
    Z  = F
  )
  
  dat <- bind_rows(dat_txblobs, dat_coblobs, dat_conear)
  
  res <- dat %>% 
    mutate(Y0 = f0_fun(X1,X2) + 
             rnorm(n(), mean=0, sd=f0_sd),
           Y1 = f0_fun(X1,X2) + tx_effect_fun(X1,X2) + 
             rnorm(n(), mean=0, sd=f0_sd),
           Y  = ifelse(Z, Y1, Y0)) %>% 
    mutate(id = 1:n(), .before=X1)
  print(paste0("SD of control outcomes: ", round(sd(res$Y0),4)))
  
  return(res)
}

if (F) {
  df <- gen_df_adv(
    nc=nc, 
    nt=nt, 
    f0_sd = f0_sd,
    tx_effect_fun = function(X1, X2) {1},
    f0_fun = function(x,y) {
      matrix(c(x,y), ncol=2) %>% 
        dmvnorm(mean = c(0.5,0.5),
                sigma = matrix(c(1,0.8,0.8,1), nrow=2))}) %>% 
    mutate(X3 = X1*X2)
  df %>% 
    ggplot(aes(X1,X2)) +
    geom_point(aes(color=Y0)) +
    facet_wrap(~Z)
}



# canonical examples ------------------------------------------------------





# old functions -----------------------------------------------------------

# simplest example in 2D, just like toy example
gen_df_adv2d <- function(nc, nt, 
                         tx_effect = 0.3,   # constant tx effect
                         sd = 0.1,
                         effect_fun = function(X1, X2) { abs(X1-X2) }) {
  
  # tx units clustered at (0,0) and (1,1)
  dat_txblobs <- tibble(
    X1 = c(rnorm(nt/2, mean=0, sd=0.2), 
           rnorm(nt/2, mean=1, sd=0.2)),
    X2 = c(rnorm(nt/2, mean=0, sd=0.2), 
           rnorm(nt/2, mean=1, sd=0.2)),
    Z  = T
  )
  
  # co units clustered at (1,0) and (0,1)
  dat_coblobs <- tibble(
    X1 = c(rnorm((nc-nt)/2, mean=0, sd=0.2), 
           rnorm((nc-nt)/2, mean=1, sd=0.2)),
    X2 = c(rnorm((nc-nt)/2, mean=1, sd=0.2), 
           rnorm((nc-nt)/2, mean=0, sd=0.2)),
    Z  = F
  )
  
  # some co units near (0,0) and (1,1)
  dat_conear <- tibble(
    X1 = c(rnorm(nt/2, mean=0, sd=0.2), 
           rnorm(nt/2, mean=1, sd=0.2)),
    X2 = c(rnorm(nt/2, mean=0, sd=0.2), 
           rnorm(nt/2, mean=1, sd=0.2)),
    Z  = F
  )
  
  dat <- bind_rows(dat_txblobs, dat_coblobs, dat_conear)
  
  res <- dat %>% 
    mutate(Y0 = effect_fun(X1,X2) + rnorm(n(), mean=0, sd=sd),
           Y1 = effect_fun(X1,X2) + tx_effect + rnorm(n(), mean=0, sd=sd),
           Y  = ifelse(Z, Y1, Y0)) %>% 
    mutate(id = 1:n(), .before=X1)
  print(paste0("SD of control outcomes: ", round(sd(res$Y0),4)))
  
  return(res)
}

# TODO: do this directly with mvnorms,
#  otherwise you get a carved-out "cross" not a carved-out "circle"
gen_df_veryadv <- function(nc, nt, 
                           tx_effect = 0.2,   # constant tx effect
                           effect_fun = function(X1, X2) { abs(X1-X2) },
                           sig = 0.2) {
  MU <- c(0.5, 0.5)
  SIG <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow=2)
  
  # tx units clustered at (0,0) and (1,1)
  dat_txblobs <- tibble(
    X1 = c(rnorm(nt/2, mean=0, sd=sig), 
           rnorm(nt/2, mean=1, sd=sig)),
    X2 = c(rnorm(nt/2, mean=0, sd=sig), 
           rnorm(nt/2, mean=1, sd=sig)),
    Z  = T
  )
  
  # co units clustered at (1,0) and (0,1)
  dat_coblobs <- tibble(
    X1 = c(rnorm((nc-nt)/2, mean=0, sd=sig), 
           rnorm((nc-nt)/2, mean=1, sd=sig)),
    X2 = c(rnorm((nc-nt)/2, mean=1, sd=sig), 
           rnorm((nc-nt)/2, mean=0, sd=sig)),
    Z  = F
  )
  
  # some co units near (0,0) and (1,1)
  dcustom <- function(x, mean, sd) {
    dnorm(x, mean=mean, sd=sd) - dnorm(x, mean=mean, sd=sd/2)/2
  }
  rcustom <- function(n, mean, sd) {
    res <- c()
    while (length(res) < n) {
      proposals <- rnorm(n, mean, 2*sd)
      accept <- runif(n) < dcustom(proposals,mean,sd) / dnorm(proposals,mean,2*sd)
      accepted <- proposals[which(accept)]
      res <- c(res, accepted)
    }
    return(res[1:n])
  }
  if (F) {
    tibble(x = seq(-5,5,0.01)) %>% 
      mutate(dens = dcustom(x, mean=0, sd=1),
             dens2 = dnorm(x, mean=0, sd=2)) %>% 
      pivot_longer(dens:dens2) %>% 
      ggplot(aes(x=x)) +
      geom_histogram(data=tibble(x=rcustom(1000, 0, 1)),
                     aes(y=..density..)) +
      geom_point(aes(y=value*2, color=name))
  }
  
  dat_conear <- tibble(
    X1 = c(rcustom(nt/2, mean=0, sd=sig), 
           rcustom(nt/2, mean=1, sd=sig)),
    X2 = c(rcustom(nt/2, mean=0, sd=sig), 
           rcustom(nt/2, mean=1, sd=sig)),
    Z  = F
  )
  
  dat <- bind_rows(dat_txblobs, dat_coblobs, dat_conear)
  
  dat %>% 
    mutate(Y0 = effect_fun(X1,X2),
           Y1 = effect_fun(X1,X2) + tx_effect,
           Y  = ifelse(Z, Y1, Y0)) %>% 
    mutate(id = 1:n(), .before=X1)
}



# simplest example in 2D, just like toy example
gen_df_advhet <- function(nc, nt, 
                          tx_effect = function(X1, X2) {(X1-0.5)^2+(X2-0.5)^2},
                          eps_sd = 0.1,
                          effect_fun = function(X1, X2) { abs(X1-X2) }) {
  
  # tx units clustered at (0,0) and (1,1)
  dat_txblobs <- tibble(
    X1 = c(rnorm(nt/2, mean=0, sd=0.2), 
           rnorm(nt/2, mean=1, sd=0.2)),
    X2 = c(rnorm(nt/2, mean=0, sd=0.2), 
           rnorm(nt/2, mean=1, sd=0.2)),
    Z  = T
  )
  
  # co units clustered at (1,0) and (0,1)
  dat_coblobs <- tibble(
    X1 = c(rnorm((nc-nt)/2, mean=0, sd=0.2), 
           rnorm((nc-nt)/2, mean=1, sd=0.2)),
    X2 = c(rnorm((nc-nt)/2, mean=1, sd=0.2), 
           rnorm((nc-nt)/2, mean=0, sd=0.2)),
    Z  = F
  )
  
  # some co units near (0,0) and (1,1)
  dat_conear <- tibble(
    X1 = c(rnorm(nt/2, mean=0, sd=0.2), 
           rnorm(nt/2, mean=1, sd=0.2)),
    X2 = c(rnorm(nt/2, mean=0, sd=0.2), 
           rnorm(nt/2, mean=1, sd=0.2)),
    Z  = F
  )
  
  dat <- bind_rows(dat_txblobs, dat_coblobs, dat_conear)
  
  res <- dat %>% 
    mutate(Y0 = effect_fun(X1,X2) + rnorm(n(), mean=0, sd=eps_sd),
           Y1 = effect_fun(X1,X2) + tx_effect(X1,X2) + rnorm(n(), mean=0, sd=eps_sd),
           Y  = ifelse(Z, Y1, Y0)) %>% 
    mutate(id = 1:n(), .before=X1)
  print(paste0("SD of control outcomes: ", round(sd(res$Y0),4)))
  
  return(res)
}




# generate co in unit square
gen_df_full <- function(nc, nt, eps_sd = 0.1,
                        tx_effect = function(X1, X2) {(X1-0.5)^2+(X2-0.5)^2},
                        effect_fun = function(X1, X2) { abs(X1-X2) }) {
  dat_co <- tibble(
    X1 = runif(nc),
    X2 = runif(nc),
    Z = F
  )
  
  dat_tx1 <- tibble(
    X1 = rnorm(nt/2, mean=0.25, sd=0.1),
    X2 = rnorm(nt/2, mean=0.25, sd=0.1),
    Z = T
  )
  dat_tx2 <- tibble(
    X1 = rnorm(nt/2, mean=0.75, sd=0.1),
    X2 = rnorm(nt/2, mean=0.75, sd=0.1),
    Z = T
  )
  
  dat <- bind_rows(dat_co, dat_tx1, dat_tx2)
  print(dat %>% mutate(eff = tx_effect(X1,X2)) %>% head())
  res <- dat %>% 
    mutate(Y0 = effect_fun(X1,X2) + rnorm(n(), mean=0, sd=eps_sd),
           Y1 = effect_fun(X1,X2) + tx_effect(X1,X2) + rnorm(n(), mean=0, sd=eps_sd),
           Y  = ifelse(Z, Y1, Y0)) %>% 
    mutate(id = 1:n(), .before=X1)
  print(paste0("SD of control outcomes: ", round(sd(res$Y0),4)))
  
  return(res)
}


# gen_df_ps <- function(n, eps_sd = 0.1,
#                       ps = function(x,y) {invlogit(2*x+2*y - 5)},
#                       tx_effect = function(X1, X2) {(X1-0.5)^2+(X2-0.5)^2},
#                       effect_fun = function(X1, X2) { abs(X1-X2) }) {
#   dat <- tibble(
#     X1 = runif(n),
#     X2 = runif(n)
#   ) %>% 
#     mutate(Z = runif(n) < ps(X1,X2))
#   
#   res <- dat %>% 
#     mutate(Y0 = effect_fun(X1,X2) + rnorm(n(), mean=0, sd=eps_sd),
#            Y1 = effect_fun(X1,X2) + tx_effect(X1,X2) + rnorm(n(), mean=0, sd=eps_sd),
#            Y  = ifelse(Z, Y1, Y0)) %>% 
#     mutate(id = 1:n(), .before=X1)
#   
#   return(res)
# }



