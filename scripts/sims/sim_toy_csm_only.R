
# Same as sim_adversarial, but it's actually the toy example used in the paper!

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

USE_PARALLEL = F
if (USE_PARALLEL) {
  library(foreach)
  library(doParallel)
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
}


# set superlearner libraries
SL.library1 <- c("SL.mean", "SL.lm", "SL.glm")

# "based on Dorie et al. ACIC 2016 competition
#   and SuperLearner documentation,
#   adjusted to reduce the computational burden"
SL.library2 <- c("SL.glm", "SL.gam", "SL.glmnet",
                 "SL.randomForest", "SL.xgboost")

# tmle defaults
SL.library3Q <- c("SL.glm", "tmle.SL.dbarts2", "SL.glmnet")   # default tmle Q.SL.library
SL.library3g <-  c("SL.glm", "tmle.SL.dbarts.k.5", "SL.gam")  # default tmle g.SL.library



# run simulations ---------------------------------------------------------

zform1 <- as.formula("Z ~ X1+X2")
zform2 <- as.formula("Z ~ X1*X2")
form1 <- as.formula("Y ~ X1+X2")
form2 <- as.formula("Y ~ X1*X2")

# repeatedly run for all combinations of pars
res <- foreach(
  i=1:40,
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

# for (i in 1:250) {
# i <- 1
  tic()

  nc <- 500
  nt <- 100
  f0_sd <- 0.5
  df <- gen_df_adv(
    nc=nc,
    nt=nt,
    f0_sd = f0_sd,
    tx_effect_fun = function(X1, X2) {3*X1+3*X2},
    # f0_fun = function(x,y) {abs(x-y)})
    f0_fun = function(x,y) {
      matrix(c(x,y), ncol=2) %>%
        dmvnorm(mean = c(0.5,0.5),
                sigma = matrix(c(1,0.8,0.8,1), nrow=2)) * 20   # multiply for more slope!
      })

  num_bins <- 6

  dist_scaling <- df %>%
    summarize(across(starts_with("X"),
                     function(x) {
                       if (is.numeric(x)) 2 * num_bins / (max(x) - min(x))
                       else 1000
                     }))



  preds_csm_avg <- get_cal_matches(
    df = df,
    metric = "maximum",
    dist_scaling = dist_scaling,
    cal_method = "fixed",
    est_method = "average",
    return = "all",
    knn = 25)

  preds_csm_scm <- get_cal_matches(
    df = df,
    metric = "maximum",
    dist_scaling = dist_scaling,
    cal_method = "fixed",
    est_method = "scm",
    return = "all",
    knn = 25)

  preds_cem <- get_cem_matches(
    df = df,
    num_bins = num_bins,
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
    nc = nc,
    nt = nt,
    f0_sd = f0_sd,
    num_bins = num_bins,
    ninf = ninf,
    ninf_cem = ninf_cem,

    true_ATT = df %>%
      filter(Z) %>%
      summarize(att = mean(Y1-Y0)) %>%
      pull(att),

    diff = get_att_diff(df),

    # bal1 = get_att_bal(df, zform1, c(0.01, 0.01)),
    # bal2 = get_att_bal(df, zform2, c(0.01, 0.01, 0.1)),
    #
    # or_lm = get_att_or_lm(df, form=form2),
    # or_bart = get_att_or_bart(df, covs=c(X1,X2)),
    # ps_lm = get_att_ps_lm(df, zform2),
    # ps_bart = get_att_ps_bart(df, covs=c(X1,X2)),
    csm_scm = get_att_point_est(preds_csm_scm),
    csm_avg = get_att_point_est(preds_csm_avg),
    # csm_scm = get_att_csm(df,
    #                       dist_scaling=dist_scaling,
    #                       est_method="scm",
    #                       rad_method = "fixed"),
    # csm_avg = get_att_csm(df,
    #                       dist_scaling=dist_scaling,
    #                       est_method="average",
    #                       rad_method = "fixed"),
    cem_scm = get_att_cem(df,
                          num_bins=num_bins,
                          est_method="scm",
                          estimand = "CEM-ATT"),
    cem_avg = get_att_cem(df,
                          num_bins=num_bins,
                          est_method="average",
                          estimand = "CEM-ATT"),
    onenn = get_att_1nn(df, dist_scaling=dist_scaling),

    # tmle1 = get_att_tmle(df,
    #                      covs=c(X1,X2),
    #                      Q.SL.library = SL.library1,
    #                      g.SL.library = SL.library1),
    # aipw1 = get_att_aipw(df,
    #                      covs=c(X1,X2),
    #                      Q.SL.library = SL.library1,
    #                      g.SL.library = SL.library1),
    #
    # tmle3 = get_att_tmle(df, covs=c(X1,X2),
    #                      Q.SL.library = SL.library3Q,
    #                      g.SL.library = SL.library3g),
    # aipw3 = get_att_aipw(df, covs=c(X1,X2),
    #                      Q.SL.library = SL.library2,
    #                      g.SL.library = SL.library2)
  )

  toc()

  # res %>%
  #   pivot_longer(-c(runid:true_ATT))

  FNAME <- "data/outputs/sim_toy_results/toy_spaceship8.csv"
  if (file.exists(FNAME)) {
    write_csv(res, FNAME, append=T)
  } else {
    write_csv(res, FNAME)
  }

  res
}

stopCluster(cl)



# analyze results ---------------------------------------------------------

# spaceship2 uses fixed calipers, not adaptive.
# spaceship3 increases slope and cuts sample size and decreases nbins
#  - cem does even better here...?
# spaceship4: more noise, heterogeneous tx effects
#  - right direction
# spaceship5: even more noise, more heterogeneous tx effects
#  - more in the right direction (csm better than cem now in terms of rmse),
#    but well-specified outcome regression lm is pretty great too
# spaceship6: double noise, less het tx effects
#  - similar? avg is still better than scm, which is annoying...
# spaceship7: more slope in f0, looks great, avg still better

# spaceship8: cut noise in half to try to get scm to do better
#  - nope, just makes onenn good. avg still better...

# spaceship7 final!


# plot sample data
set.seed(1)
df <- gen_df_adv(
  nc=500,
  nt=100,
  f0_sd = 0.1,
  tx_effect_fun = function(X1, X2) {0.2},
  # f0_fun = function(x,y) {abs(x-y)})
  f0_fun = function(x,y) {
    matrix(c(x,y), ncol=2) %>%
      dmvnorm(mean = c(0.5,0.5),
              sigma = matrix(c(1,0.8,0.8,1), nrow=2)) * 10   # multiply for more slope!
  })
vlines <- seq(min(df$X1),max(df$X1),length.out=9)
hlines <- seq(min(df$X2),max(df$X2),length.out=9)
df %>%
  ggplot(aes(x=X1, y=X2, color=Y)) +
  geom_point(aes(pch=Z)) +
  scale_color_continuous(low="blue", high="orange") +
  geom_vline(xintercept=vlines, lty="dotted") +
  geom_hline(yintercept=hlines, lty="dotted") +
  facet_wrap(~Z)

d1 <- (vlines[2]-vlines[1])/2
d2 <- (hlines[2]-hlines[1])/2
df %>%
  ggplot(aes(x=X1, y=X2, color=Y)) +
  geom_rect(data=. %>% filter(Z) %>% slice(1:5),
            aes(xmin=X1-d1, xmax=X1+d1, ymin=X2-d2, ymax=X2+d2),
            fill = NA, color="black", alpha=0.1) +
  geom_point(aes(pch=Z)) +
  scale_color_continuous(low="blue", high="orange") +
  geom_vline(xintercept=vlines, lty="dotted") +
  geom_hline(yintercept=hlines, lty="dotted")

if (F) {
  # for seed(1): cem does better! it drops a bunch though.
  df %>%
    filter(Z) %>%
    summarize(ATT = mean(Y1-Y0))
  csm_scm <- get_cal_matches(
    df = df,
    metric = "maximum",
    cal_method = "adaptive",
    est_method = "scm",
    dist_scaling = df %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) 8 / (max(x) - min(x))
                         else 1000
                       })),
    return = "all",
    knn = 25)
  csm_scm %>%
    count(subclass) %>%
    ggplot(aes(x=n)) +
    geom_histogram(binwidth=2, color="black")

  csm_avg = get_att_csm(df, num_bins=8, est_method="average")
  cem_scm = get_att_cem(df, num_bins=8, est_method="scm")
  cem_avg = get_att_cem(df, num_bins=8, est_method="average")

  csm_scm
  csm_avg
  cem_scm
  cem_avg
}




total_res <- read_csv("data/outputs/sim_toy_results/toy_spaceship8.csv")

res_df <- total_res %>%
  pivot_longer(-(runid:true_ATT)) %>%
  group_by(nc,nt,f0_sd,name) %>%
  summarize(rmse = sqrt(mean((value-true_ATT)^2)),
            bias = mean(value-true_ATT),
            sd   = sd(value)) %>%
  rename(method = name)

res_df %>%
  pivot_longer(c(rmse, bias, sd)) %>%
  ggplot(aes(x=method)) +
  geom_col(aes(y=value)) +
  facet_wrap(~name, scales="free") +
  coord_flip()

# prove identity
total_res %>%
  pivot_longer(-(runid:true_ATT)) %>%
  group_by(nc,nt,f0_sd,name) %>%
  summarize(rmse = sqrt(mean((value-true_ATT)^2)),
            bias = mean(value-true_ATT),
            var_est = mean((value-mean(value))^2),
            var_true = mean((true_ATT-mean(true_ATT))^2),
            cov = mean((value-mean(value))*(true_ATT-mean(true_ATT)))) %>%
  # mutate(rmse2 = sqrt(var_est + var_true + bias^2 - 2*cov),
  #        correct = all.equal(rmse,rmse2)) %>%
  rename(method = name) %>%
  pivot_longer(c(rmse, bias, var_est, var_true, cov)) %>%
  ggplot(aes(x=method)) +
  geom_col(aes(y=value)) +
  facet_wrap(~name, scales="free") +
  coord_flip()

# results:
#  - toy_constant-tx-effect:
#     - aipw1/bal2/csm_avg best...?
#     - 1nn lowest bias but very variable
#     - in general, does seem like lowering variance matters quite a bit...
#  - toy_constant-tx-effect2: doubled dnorm f0 function, for more slope there & cut noise
#     - 1nn and cem_scm suddenly very good, csm better than bal2/aipw1 but not great...
#     --> I think 1nn does very well with low noise.
#  - toy_constant-tx-effect3: triple dnorm, noise back up to 0.1
#     - still, 1nn and cem-scm doing very well in terms of bias for whatever reason...
#     - but csm better than everything else! moving in the right direction...?
#  - toy_constant-tx-effects4: quadruple dnorm, noise up to 0.2, smaller cells, more tx units
#     - pretty similar to 3...

# Q: HOW in the world does CEM outperform CSM in terms of bias?

#  - toy_constant-tx-effects5: same as 4, but fixed bug where
#    csm is calipering on X3 as well (annoying that this changes things),
#    also added X3 to tmle estimator
#     - results: csm better now in terms of rmse, tmle1 not horrible anymore.
#     - challenges: avg is better than scm in terms of bias??
#                   cem is better than csm in terms of bias??
#  - toy_constant-tx-effects6: reduce noise again?
#     - doesn't seem to change results








