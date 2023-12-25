
# functions for simulating data

require(aciccomp2016)


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

# generate sample dataset from Hainmueller (2012), exactly
#  - note: generates a big population and samples nt/nco units
# paper settings:
#  - n (total number of units): 300, 600, 1500
#  - r (ratio of co/tx units): 1,2,5
#  - sigma_e: n30 = low overlap, n100 = high overlap, chi5 = weird overlap
#  - outcome: "linear", "nl1", "nl2"
#  - sigma_y: 1
gen_df_hain <- function(nt = 50,
                        nc = 250,
                        sigma_e = c("chi5", "n30", "n100"), 
                        outcome = c("linear", "nl1", "nl2"),
                        sigma_y = 1,
                        ATE = 0) {
  sigma_e <- match.arg(sigma_e)
  outcome <- match.arg(outcome)
  NUMSAMP <- max(nt, nc)*10
  
  if (sigma_e == "n30") {
    eps_e <- rnorm(NUMSAMP, sd=sqrt(30))
  } else if (sigma_e == "n100"){
    eps_e <- rnorm(NUMSAMP, sd=sqrt(100))
  } else if (sigma_e == "chi5") {
    eps_e <- rchisq(NUMSAMP, df=5)
    eps_e <- (eps_e - mean(eps_e)) / sd(eps_e)
    eps_e <- eps_e * sqrt(67.6) + 0.5
  }
  
  df <- as_tibble(rmvnorm(NUMSAMP, 
                          mean = c(0,0,0), 
                          sigma = matrix(c(2,1,-1,1,1,-0.5,-1,-0.5,1), ncol=3))) %>%
    mutate(V4 = runif(NUMSAMP, min=-3, max=3),
           V5 = rchisq(NUMSAMP, df=1),
           V6 = sample(c(T,F), NUMSAMP, replace=T, prob=c(0.5,0.5)),
           Z  = (V1 + 2*V2 - 2*V3 + V4 - 0.5*V5 + V6 + eps_e) > 0,
           id = 1:n()) %>%
    relocate(id)
  
  if (outcome == "linear") {
    df <- df %>%
      mutate(Y0 = V1 + V2 + V3 - V4 + V5 + V6 + rnorm(NUMSAMP, sd=sigma_y))
  } else if (outcome == "nl1") {
    df <- df %>%
      mutate(Y0 = V1 + V2 + 0.2*V3*V4 - sqrt(V5) + rnorm(NUMSAMP, sd=sigma_y))
  } else if (outcome == "nl2") {
    df <- df %>%
      mutate(Y0 = (V1 + V2 + V5)^2 + rnorm(NUMSAMP, sd=sigma_y))
  } else if (outcome == "nl3") {   # not in paper
    df <- df %>%
      mutate(Y0 = V1*V2 + V2*V5 + V1*V5 + 3*cos(V3) + 2*sin(V2) + 
               V3*V4 + 0.5*V4^2 + rnorm(NUMSAMP, sd=sigma_y))
  }
  df <- df %>%
    mutate(Y1 = Y0 + ATE,
           Y  = ifelse(Z, Y1, Y0))
  
  # get correct number of observations
  if (nt > sum(df$Z) | nc > nrow(df)-sum(df$Z)) {
    warning("Insufficient population pool: increase NUMSAMP")
  }
  df_tx <- df %>%
    filter(Z) %>%
    slice_sample(n=nt)
  df_co <- df %>%
    filter(!Z) %>%
    slice_sample(n=nc)
  
  # bind rows, tx then co
  df_final <- df_tx %>%
    bind_rows(df_co)
  names(df_final) <- names(df_final) %>% 
    str_replace("V", "X")
  
  return(df_final)
}

gen_df_acic <- function(model.trt="step", 
                        root.trt=0.35, 
                        overlap.trt="full",
                        model.rsp="linear", 
                        alignment=0.75, 
                        te.hetero="high",
                        random.seed=1,
                        n=1000, p=10) {
  # idea: input_2016 is a 4802 x 58 df of COVARIATES
  #  - dgp_2016() outputs a dataset of POs for input_2016
  #  - see parameters_2016 for example inputs
  
  # rows:
  #  - model.trt: tx assignment model (linear, polynomial, step)
  #  - root.trt: percentage of tx units (any #; 0.35, 0.65)
  #  - overlap.trt: overlap between tx/co (one-term, full)
  #  - model.rsp: outcome model (linear, exponential, step)
  #  - alignment: prop of tx assignment covs used in response fxn (any #; 0, 0.25, 0.75)
  #  - te.hetero: tx-effect heterogeneity (none, med, high)
  
  # NOTE: for some god-forsaken reason, subsetting rows breaks things
  #  - i.e., it makes the browser() functionality start running...
  #  - the browser() call is from newDiscreteBaseFunction(),
  #    where numLevels <= 1...?
  
  set.seed(random.seed)

  # only subset numeric columns
  numeric_cols <- input_2016 %>% 
    summarize(across(everything(), ~length(unique(.)))) %>% 
    pivot_longer(everything()) %>% 
    filter(value >= 20) %>% 
    pull(name)
  base_df <- input_2016[,numeric_cols]
  
  if (n > nrow(base_df)) {
    error(glue("n too large; needs to be <= {nrow(base_df)}"))
  }
  if (p > ncol(base_df)) {
    error(glue("p too large; needs to be <= {ncol(base_df)}"))
  }
  base_df <- base_df[,1:p] %>%
    sample_n(n, replace=F)
  
  # # NECESSARY: ensure that there aren't any factor levels without entries
  # base_df$x_2 <- as.factor(as.character(base_df$x_2))
  
  # base_df %>%
  #   summarize(across(everything(), ~length(unique(.))))
  # browser()
  df <- 
    dgp_2016(base_df,
             parameters = list(
               model.trt = model.trt,
               root.trt = root.trt,
               overlap.trt = overlap.trt,
               model.rsp = model.rsp,
               alignment = alignment,
               te.hetero = te.hetero),
             random.seed = random.seed)
  
  df <- bind_cols(base_df, df) %>%
    rename_with(.cols = everything(), toupper) %>%
    rename_with(.cols = everything(), ~str_remove(., "_")) %>%
    rename_with(.cols = everything(), ~str_remove(., "\\.")) %>%
    mutate(id = 1:n(), .before=X1,
           Z = Z==1) %>%
    as_tibble()
  
  return(df)
}

gen_df_kang <- function(n=1000) {
  as_tibble(rmvnorm(n, mean=rep(0,4), sigma=diag(4))) %>% 
    mutate(Y = 210 + 27.4*V1 + 13.7*(V2+V3+V4) + rnorm(n),
           e = expit(-V1 + 0.5*V2 - 0.25*V3 - 0.1*V4)) %>% 
    mutate(X1 = exp(V1/2),
           X2 = V2/(1+exp(V1)) + 10,
           X3 = (V1*V3/25 + 0.6)^3,
           X4 = (V2 + V4 + 20)^2) %>% 
    rowwise() %>% 
    mutate(Z = sample(c(T,F), size=1, replace=T, 
                      prob=c(e,1-e))) %>% 
    ungroup() %>% 
    mutate(id = 1:n, .before=V1)
}


if (F) {
  gen_df_acic(n=1000, p=10)
  
  input_2016[,1:10] %>% 
    sample_n(1000,replace=F) %>% 
    dgp_2016(
      parameters = list(
        model.trt="linear", 
        root.trt=0.35, 
        overlap.trt="full", 
        model.rsp="linear", 
        alignment=0.75, 
        te.hetero="none"), 
        random.seed=1)
}



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


### Generate Otsu and Rai DGP to debug the bootstrap
gen_df_otsu <- function(N = 1000,K = 2){
  # output: "id" "X1" "X2" "Z"  "Y0" "Y1" "Y" 
  # Function m(z) as defined
  # K=2; N = 1000
  m <- function(z) {
    0.4 + 0.25 * sin(8 * z - 5) + 0.4 * exp(-16 * (4 * z - 2.5)^2)
  }
  
  # Generate data
  xi <- runif(N, 0, 1)
  zeta <- matrix(rnorm(N * K), ncol = K)
  X <- xi * apply(abs(zeta), 1, function(x) x / sqrt(sum(x^2)))
  norm_X <- sqrt(colSums(X^2))
  
  # Parameters for P(X)
  gamma1 <- 0.15
  gamma2 <- 0.7
  P_X <- gamma1 + gamma2 * norm_X
  
  # Treatment assignment
  vi <- runif(N)
  Z <- as.numeric(P_X >= vi)
  
  # Potential outcomes
  epsilon <- rnorm(N, 0, 0.2)
  Y0 <- m(norm_X) + epsilon
  Y1 <- Y0 + 0  # tau is set to 0
  Y <- (1 - Z) * Y0 + Z * Y1
  
  # Create dataframe
  df <- data.frame(id = 1:N, 
                   X1 = X[1, ], 
                   X2 = X[2, ], 
                   Z = Z, 
                   Y0 = Y0, 
                   Y1 = Y1, 
                   Y = Y)
  return(df)
}

generate_one_otsu <- function(){
  gen_df_otsu(N = 1000,K = 2)
}

generate_one_toy <- function(){
  gen_df_adv(
    nc=500, 
    nt=100, 
    f0_sd = 0.5,
    tx_effect_fun = function(X1, X2) {3*X1+3*X2},
    f0_fun = function(x,y) {
      matrix(c(x,y), ncol=2) %>%
        dmvnorm(mean = c(0.5,0.5),
                sigma = matrix(c(1,0.8,0.8,1), nrow=2)) * 20   # multiply for more slope!
    })
}

get_df_scaling_from_dgp_name <- 
  function(dgp_name,
           kang_true = F){
  if (dgp_name == "toy"){
    df_dgp <- generate_one_toy()
    dist_scaling <- df_dgp %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) 6 / (max(x) - min(x))
                         else 1000
                       }))
  }else if(dgp_name=="kang"){
    df_dgp <- gen_df_kang(n = 1000)
    if (kang_true){ # Replace X1 to X4 by V1 to V4
      df_dgp[, paste0("X",1:4)] <-
        df_dgp %>% select(starts_with("V"))
    }
    dist_scaling <- df_dgp %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) 5 / (max(x) - min(x))
                         else 1000
                       }))
  }else if (dgp_name == "otsu"){
    df_dgp <- generate_one_otsu()
    dist_scaling <- df_dgp %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) 6 / (max(x) - min(x))
                         else 1000
                       }))
  } else {
    stop("dgp_name must be toy or kang or otsu")
  }
  return(list(df_dgp=df_dgp,
              dist_scaling=dist_scaling))
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



