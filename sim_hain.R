
# Hainmueller simulation

require(tidyverse)
require(mvtnorm)

require(optweight)
require(dbarts)
require(tmle)
require(AIPW)

require(tictoc)

source("R/distance.R")
source("R/sc.R")
source("R/matching.R")
source("R/estimate.R")
source("R/sim_data.R")
source("R/wrappers.R")
source("R/utils.R")


# set superlearner libraries
SL.library1 <- c("SL.mean", "SL.lm", "SL.glm")

# "based on Dorie et al. ACIC 2016 competition
#   and SuperLearner documentation,
#   adjusted to reduce the computational burden"
SL.library2 <- c("SL.glm", #"SL.gam", 
                 "SL.glmnet",
                 # "SL.gbm", "SL.bartMachine",   # too slow
                 "SL.randomForest", "SL.xgboost")             # based on aipw recommendations

# tmle defaults
SL.library3Q <- c("SL.glm", "tmle.SL.dbarts2", "SL.glmnet")   # default tmle Q.SL.library
SL.library3g <-  c("SL.glm", "tmle.SL.dbarts.k.5", "SL.gam")  # default tmle g.SL.library



# run simulations ---------------------------------------------------------

zform1 <- as.formula("Z ~ X1+X2+X3+X4+X5+X6")
zform2 <- as.formula("Z ~ (X1+X2+X3+X4+X5+X6)^2")
form1 <- as.formula("Y ~ X1+X2+X3+X4+X5+X6")
form2 <- as.formula("Y ~ (X1+X2+X3+X4+X5+X6)^2")

# repeatedly run for all combinations of pars
for (i in 1:39) {
  nc <- 250
  nt <- 50
  df <- gen_df_hain(
    nc=nc, 
    nt=nt, 
    sigma_e = "n100",   # high overlap condition
    outcome = "nl2",
    sigma_y = 1,
    ATE = 0)
  
  nbins <- 5
  dist_scaling <- df %>%
    summarize(across(starts_with("X"),
                     function(x) {
                       if (is.numeric(x)) nbins / (max(x) - min(x))
                       else 1000
                     }))
  
  tic()
  
  # for hain simulation: drop units that don't get a within-caliper matches
  preds_csm <- get_cal_matches(
    df = df,
    metric = "maximum",
    dist_scaling = dist_scaling,
    cal_method = "fixed",
    est_method = "average",
    return = "all",
    knn = 25)
  preds_cem <- get_cem_matches(
    df = df,
    num_bins = nbins,
    est_method = "average",
    return = "all")
  
  ninf <- length(attr(preds_csm, "unmatched_units"))
  ninf_cem <- nt - length(attr(preds_cem, "feasible_units"))
  
  # for all methods: filter out unmatched csm units
  df <- df %>% 
    filter(!id %in% attr(preds_csm, "unmatched_units"))
  
  res <- tibble(
    runid = i,
    ninf = ninf,
    ninf_cem = ninf_cem,
    true_ATT = 0,
    
    diff = get_att_diff(df),
    
    bal1 = get_att_bal(df, zform1, rep(0.01, 6)),
    bal2 = get_att_bal(df, zform2, rep(0.01, 21)),
    
    or_lm = get_att_or_lm(df, form=form2),
    or_bart = get_att_or_bart(df, covs=c(X1,X2,X3,X4,X5,X6)),
    ps_lm = get_att_ps_lm(df, zform2),
    ps_bart = get_att_ps_bart(df, covs=c(X1,X2,X3,X4,X5,X6)),
    
    csm_scm = get_att_csm(df, dist_scaling=dist_scaling, est_method="scm",
                          cal_method = "fixed"),
    csm_avg = get_att_csm(df, dist_scaling=dist_scaling, est_method="average",
                          cal_method = "fixed"),
    cem_scm = get_att_cem(df, num_bins=nbins, est_method="scm",
                          estimand = "CEM-ATT"),
    cem_avg = get_att_cem(df, num_bins=nbins, est_method="average",
                          estimand = "CEM-ATT"),
    onenn = get_att_1nn(df, dist_scaling=dist_scaling),
    
    tmle1 = get_att_tmle(df, covs=c(X1,X2,X3,X4,X5,X6),
                         Q.SL.library = SL.library1,
                         g.SL.library = SL.library1),
    tmle2 = get_att_tmle(df, covs=c(X1,X2,X3,X4,X5,X6),
                         Q.SL.library = SL.library3Q,
                         g.SL.library = SL.library3g),
    
    aipw1 = get_att_aipw(df, covs=c(X1,X2,X3,X4,X5,X6),
                         Q.SL.library = SL.library1,
                         g.SL.library = SL.library1),
    aipw2 = get_att_aipw(df, covs=c(X1,X2,X3,X4,X5,X6),
                         Q.SL.library = SL.library2,
                         g.SL.library = SL.library2)
  )
  
  toc()
  
  res %>% 
    pivot_longer(-c(runid:true_ATT))
  
  FNAME <- "sim_canonical_results/hain_final.csv"
  if (file.exists(FNAME)) {
    write_csv(res, FNAME, append=T)
  } else {
    write_csv(res, FNAME)
  }
}


# analyze results ---------------------------------------------------------

# hain_testing: try to replicate Hain
#  - I think Table 1 in Hain uses MSE * 100?
# hain_testing2: try to replicate my old results
#  - not at all...? old results are far better...???
#  - ISSUE: adaptive caliper extrapolation is v bad for the estimate

# OPEN Q: why is this simulation showing CSM so much worse than 
#  my sim_hain analyses in scmatch1???

# TODO:
#  - I think that Table 1 in Hain uses MSE * 100? Seems similar for 1nn
#  - Goal: show that csm isn't too bad... it's pretty bad...

total_res <- read_csv("sim_canonical_results/hain_final.csv")

total_res %>% 
  pivot_longer(-(runid:true_ATT)) %>%
  group_by(name) %>%
  summarize(# mse = mean((value-true_ATT)^2),
            rmse = sqrt(mean((value-true_ATT)^2)),
            bias = mean(value-true_ATT),
            sd   = sd(value)) %>% 
  rename(method = name) %>% 
  pivot_longer(c(# mse, 
                 rmse, bias, sd)) %>%
  ggplot(aes(x=method)) +
  geom_col(aes(y=value)) +
  facet_wrap(~name, scales="free") +
  coord_flip()

total_res %>% 
  pivot_longer(-(runid:true_ATT)) %>%
  group_by(name) %>%
  summarize(# mse = mean((value-true_ATT)^2),
    rmse = sqrt(mean((value-true_ATT)^2)),
    bias = mean(value-true_ATT),
    sd   = sd(value)) %>% 
  rename(method = name) %>% 
  filter(method %in% c("tmle1", "tmle2", "aipw1", "aipw2", "onenn", "csm_scm", "cem_avg", "bal1")) %>% 
  pivot_longer(c(# mse, 
    rmse, bias, sd)) %>%
  ggplot(aes(x=method)) +
  geom_col(aes(y=value)) +
  facet_wrap(~name, scales="free") +
  coord_flip()




