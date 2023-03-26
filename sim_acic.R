
# ACIC simulation

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

# parallelize
if (T) {
  library(foreach)
  library(doParallel)
  cores=detectCores()
  cl <- makeCluster(cores[1]-1, outfile="") #not to overload your computer
  registerDoParallel(cl)
  # registerDoSEQ()   # sequential, for easier debugging
}

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

tic()


# run simulations ---------------------------------------------------------

# repeatedly run for all combinations of pars
# for (i in 1:100) {
res <- foreach(
  i=sample(1:100000, 100),    # random seeds
  .packages = c("tidyverse", 
                "mvtnorm", "optweight", "dbarts", "tmle", "AIPW", "tictoc",
                "aciccomp2016"),
  .combine=rbind) %dopar% {
    source("R/distance.R")
    source("R/sc.R")
    source("R/matching.R")
    source("R/estimate.R")
    source("R/sim_data.R")
    source("R/wrappers.R")
    source("R/utils.R")
    
  n <- 1000
  p <- 10
  model.trt = "step"
  model.rsp = "step"
  
  df <- gen_df_acic(
    model.trt=model.trt, 
    root.trt=0.35, 
    overlap.trt="full",
    model.rsp=model.rsp, 
    alignment=0.75, 
    te.hetero="high",
    random.seed=i,        # NOTE: set random seed!
    n=n, 
    p=p
  )
  
  covs <- df %>% 
    select(starts_with("X")) %>% 
    colnames()
  zform1 <- as.formula(paste0("Z ~ ", paste0(covs, collapse="+")))
  zform2 <- as.formula(paste0("Z ~ (", paste0(covs, collapse="+"), ")^2"))
  form1 <- as.formula(paste0("Y ~ ", paste0(covs, collapse="+")))
  form2 <- as.formula(paste0("Y ~ (", paste0(covs, collapse="+"), ")^2"))
  
  nbins <- 5
  dist_scaling <- df %>%
    summarize(across(starts_with("X"),
                     function(x) {
                       if (is.numeric(x)) nbins / (max(x) - min(x))
                       else 1000
                     }))
  
  tic()
  
  # for acic simulation: drop units that don't get a within-caliper matches
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
  
  # record number of infeasible units
  ninf <- length(attr(preds_csm, "unmatched_units"))
  ninf_cem <- sum(df$Z) - length(attr(preds_cem, "feasible_units"))
  
  # for all methods: filter out unmatched csm units
  df <- df %>% 
    filter(!id %in% attr(preds_csm, "unmatched_units"))
  
  res <- tibble(
    runid = i,
    n = n,
    p = p,
    model.trt = model.trt,
    model.rsp = model.rsp,
    ninf = ninf,
    ninf_cem = ninf_cem,
    true_ATT = df %>% 
      filter(Z) %>%
      summarize(att = mean(Y1-Y0)) %>%
      pull(att),
    
    diff = get_att_diff(df),
    
    bal1 = get_att_bal(df, zform1, rep(0.01, length(covs))),
    bal2 = get_att_bal(df, zform2, 
                       rep(0.1, length(covs) + choose(length(covs), 2))),
    
    or_lm = get_att_or_lm(df, form=form2),
    or_bart = get_att_or_bart(df, covs=covs),
    ps_lm = get_att_ps_lm(df, zform2),
    ps_bart = get_att_ps_bart(df, covs=covs),
    
    csm_scm = get_att_csm(df, dist_scaling=dist_scaling, est_method="scm",
                          cal_method = "fixed"),
    csm_avg = get_att_csm(df, dist_scaling=dist_scaling, est_method="average",
                          cal_method = "fixed"),
    cem_scm = get_att_cem(df, num_bins=nbins, est_method="scm",
                          estimand = "CEM-ATT"),
    cem_avg = get_att_cem(df, num_bins=nbins, est_method="average",
                          estimand = "CEM-ATT"),
    onenn = get_att_1nn(df, dist_scaling=dist_scaling),
    
    tmle1 = get_att_tmle(df, covs=covs,
                         Q.SL.library = SL.library1,
                         g.SL.library = SL.library1),
    aipw1 = get_att_aipw(df, covs=covs,
                         Q.SL.library = SL.library1,
                         g.SL.library = SL.library1),

    tmle2 = get_att_tmle(df, covs=covs,
                         Q.SL.library = SL.library3Q,
                         g.SL.library = SL.library3g),
    aipw2 = get_att_aipw(df, covs=covs,
                         Q.SL.library = SL.library2,
                         g.SL.library = SL.library2)
  )
  
  toc()
  
  res %>% 
    pivot_longer(-c(runid:true_ATT))
  
  FNAME <- "sim_canonical_results/acic_spaceship.csv"
  if (file.exists(FNAME)) {
    write_csv(res, FNAME, append=T)
  } else {
    write_csv(res, FNAME)
  }
  }

toc()


# analyze results ---------------------------------------------------------


# acic_test: model.trt = "polynomial", model.rsp = "linear
#  - cem_scm okay (better than tmle1/aipw1/bal1), 
#     but slightly worse than the lm models...
# acic_test2: step, step
#  - now, cem_scm slightly better than lm models


total_res <- read_csv("sim_canonical_results/acic_spaceship.csv")

total_res %>% 
  ggplot() +
  geom_density(aes(x=ninf)) +
  geom_density(aes(x=ninf_cem), color="red")

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




