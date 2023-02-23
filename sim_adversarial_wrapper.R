
# Idea: simplest possible simulation study where:
#  1) balancing performs poorly (if not balanced on enough things)
#  2) CSM performs on part with BART

# TODO: implement more methods
#  - adaptive hyperboxes
#  - doubly robust something

# TODO: another simulation study where nothing works, 
#  but seeing which methods flag this (and what estimates get produced)

# TODO: a more complex version of this study
#  - heterogeneous treatment effects?
#  - more covariates that interact like this?
#  - open Q: what would the goal of this more complex simulation be?
#     we've already shown (in this simulation) that there are cases 
#     where CSM is nice...

require(tidyverse)
require(mvtnorm)
require(dbarts)
require(optweight)
require(tictoc)
source("sim_functions.R")



# functions for generating data -------------------------------------------

# simplest example in 2D, just like toy example
#  - constant treatment effect
gen_df_adv2d <- function(nc, nt, sig=0.2) {
  MU <- c(0.5, 0.5)
  SIG <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow=2)
  
  # tx units clustered at (0,0) and (1,1)
  dat_txblobs <- tibble(
    X1 = c(rnorm(nt, mean=0, sd=sig), rnorm(nt, mean=1, sd=sig)),
    X2 = c(rnorm(nt, mean=0, sd=sig), rnorm(nt, mean=1, sd=sig)),
    Z  = T
  )
  
  # co units clustered at (1,0) and (0,1)
  dat_coblobs <- tibble(
    X1 = c(rnorm(nc, mean=0, sd=sig), rnorm(nc, mean=1, sd=sig)),
    X2 = c(rnorm(nc, mean=1, sd=sig), rnorm(nc, mean=0, sd=sig)),
    Z  = F
  )
  
  # some co units near (0,0) and (1,1)
  dat_conear <- tibble(
    X1 = c(rnorm(nt, mean=0, sd=sig), rnorm(nt, mean=1, sd=sig)),
    X2 = c(rnorm(nt, mean=0, sd=sig), rnorm(nt, mean=1, sd=sig)),
    Z  = F
  )
  
  dat <- bind_rows(dat_txblobs, dat_coblobs, dat_conear)
  
  # function for control POs
  effect_fun <- function(X1, X2) { abs(X1-X2) }
  # constant tx effect
  tx_effect <- 0.2
  
  dat %>% 
    mutate(Y0 = effect_fun(X1,X2),
           Y1 = effect_fun(X1,X2) + tx_effect,
           Y  = ifelse(Z, Y1, Y0)) %>% 
    mutate(id = 1:n(), .before=X1)
  # rowwise() %>% 
  # mutate(Y0 = dmvnorm(x = c(X1,X2),
  #                     mean = MU,
  #                     sigma = SIG)) %>% 
  # ungroup()
}


# wrapper functions for getting att estimate ------------------------------

get_att_bal <- function(d, 
                        form,
                        tols) {
  m_bal <- optweight(form,
                     data = d,
                     tols = c(0.01, 0.01),
                     estimand = "ATT")
  
  # output ATT estimate
  d %>% 
    mutate(wt = m_bal$weights) %>% 
    group_by(Z) %>% 
    summarize(Y = mean(Y*wt)) %>% 
    pivot_wider(names_from=Z, names_prefix="Y", values_from=Y) %>% 
    mutate(ATT = YTRUE - YFALSE) %>% 
    pull(ATT)
}

get_att_bart <- function(d,
                         covs,
                         outcome) {
  
  m_bart <- bart(x.train = d %>% select(Z,{{covs}}), 
                 y.train = d %>% pull({{outcome}}),
                 x.test  = d %>% 
                   select(Z,{{covs}}) %>% 
                   filter(Z==1) %>% 
                   mutate(Z=0) )
  
  # output ATT estimate
  d %>% 
    filter(Z==1) %>% 
    mutate(Y0hat = colMeans(m_bart$yhat.test)) %>% 
    summarize(ATThat = mean(Y-Y0hat)) %>% 
    pull(ATThat)
}

get_att_lm <- function(d,
                       form) {
  m_lm <- lm(form, data=d)
  
  d %>% 
    filter(Z==T) %>% 
    mutate(Y0hat = predict(m_lm, 
                           newdata = d %>% 
                             filter(Z==T) %>% 
                             mutate(Z=F))) %>% 
    summarize(ATThat = mean(Y-Y0hat)) %>% 
    pull(ATThat)
}

get_att_csm <- function(d,
                        num_bins) {
  preds_csm <- get_cal_matches(
    df = d,
    metric = "maximum",
    cal_method = "cem",
    est_method = "scm",
    return = "sc_units",
    num_bins = num_bins,
    wider = F,
    cem_tx_units = d %>% 
      filter(Z==1) %>% 
      pull(id)
  )
  
  # get ATT estimate:
  preds_csm %>% 
    group_by(Z) %>% 
    summarize(Y = mean(Y*weights)) %>% 
    pivot_wider(names_from=Z, names_prefix="Y", values_from=Y) %>% 
    mutate(ATT = YTRUE - YFALSE) %>% 
    pull(ATT)
}

get_att_cem <- function(d,
                        num_bins) {
  preds_cem <- get_cem_matches(
    df = d,
    num_bins = 5,
    method = "average",
    return = "sc_units")
  
  # get ATT estimate:
  preds_cem %>% 
    group_by(Z) %>% 
    summarize(Y = mean(Y*weights)) %>% 
    pivot_wider(names_from=Z, names_prefix="Y", values_from=Y) %>% 
    mutate(ATT = YTRUE - YFALSE) %>% 
    pull(ATT)
}


# test these
if (F) {
  df <- gen_df_adv2d(nc=1000, nt=10, sig=0.2)
  get_att_bal(df, 
              form=as.formula('Z ~ X1+X2'), 
              tols=c(0.01, 0.01))
  get_att_bart(df, covs=c(X1, X2), outcome=Y)
  get_att_lm(df, as.formula('Y~Z+X1*X2'))
  get_att_csm(df, num_bins=5)
  get_att_cem(df, num_bins=5)
}



# run simulations ---------------------------------------------------------

pars <- expand_grid(
  nc = 1000,
  nt = 10,
  sig = 0.2
)

# repeatedly run for all combinations of pars
tic()
total_res <- map_dfr(1:250,
                     function(x) {
                       pars %>% 
                         mutate(runid = x, .before=nc) %>% 
                         rowwise() %>% 
                         mutate(df = list(gen_df_adv2d(nc, nt, sig))) %>% 
                         mutate(bal = get_att_bal(df, as.formula('Z ~ X1+X2'), c(0.01, 0.01)),
                                bart = get_att_bart(df, covs=c(X1, X2), outcome=Y),
                                lm  = get_att_lm(df, as.formula('Y~Z+X1*X2')),
                                csm = get_att_csm(df, num_bins=5),
                                cem = get_att_cem(df, num_bins=5)) %>% 
                         ungroup() %>% 
                         select(-df)
                     })
toc()

write_csv(total_res, "sim_adversarial_results/test2.csv")



# analyze results ---------------------------------------------------------

# plot sample data
gen_df_adv2d(nc=1000, nt=10, sig=0.2) %>% 
  ggplot(aes(x=X1, y=X2, color=Y)) +
  geom_point(aes(pch=Z)) +
  scale_color_continuous(low="blue", high="orange")


total_res <- read_csv("test2.csv")
# show results
total_res %>% 
  pivot_longer(bal:cem) %>% 
  group_by(name) %>% 
  summarize(rmse = sqrt(mean((value-0.2)^2)),
            bias = mean(value-0.2),
            sd   = sd(value))




