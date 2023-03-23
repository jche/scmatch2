
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


# set superlearner libraries
SL.library1 <- c("SL.mean", "SL.lm", "SL.glm")

# "based on Dorie et al. ACIC 2016 competition
#   and SuperLearner documentation,
#   adjusted to reduce the computational burden"
SL.library2 <- c("SL.glm", "SL.gam", "SL.glmnet",
                 # "SL.gbm", "SL.bartMachine",   # too slow
                 "SL.randomForest", "SL.xgboost")             # based on aipw recommendations

# tmle defaults
SL.library3Q <- c("SL.glm", "tmle.SL.dbarts2", "SL.glmnet")   # default tmle Q.SL.library
SL.library3g <-  c("SL.glm", "tmle.SL.dbarts.k.5", "SL.gam")  # default tmle g.SL.library



# run simulations ---------------------------------------------------------

zform1 <- as.formula("Z ~ X1+X2")
zform2 <- as.formula("Z ~ X1*X2")
form1 <- as.formula("Y ~ X1+X2")
form2 <- as.formula("Y ~ X1*X2")

# repeatedly run for all combinations of pars
tic()
for (i in 1:250) {
  nc <- 1000
  nt <- 100
  f0_sd <- 0.1
  df <- gen_df_adv(
    nc=nc, 
    nt=nt, 
    f0_sd = f0_sd,
    tx_effect_fun = function(X1, X2) {0.2},
    f0_fun = function(x,y) {
      matrix(c(x,y), ncol=2) %>% 
        dmvnorm(mean = c(0.5,0.5),
                sigma = matrix(c(1,0.8,0.8,1), nrow=2))}) %>% 
    mutate(X3 = X1*X2)
  
  res <- tibble(
    runid = i,
    nc = nc,
    nt = nt,
    f0_sd = f0_sd,

    true_ATT = df %>%
      filter(Z) %>%
      summarize(att = mean(Y1-Y0)) %>%
      pull(att),

    bal1 = get_att_bal(df, zform1, c(0.01, 0.01)),
    bal2 = get_att_bal(df, zform2, c(0.01, 0.01, 0.1)),

    or_lm = get_att_or_lm(df, form=form2),
    or_bart = get_att_or_bart(df, covs=c(X1,X2)),
    ps_lm = get_att_ps_lm(df, zform2),
    ps_bart = get_att_ps_bart(df, covs=c(X1,X2)),
    
    csm_scm = get_att_csm(df, num_bins=5, est_method="scm"),
    csm_avg = get_att_csm(df, num_bins=5, est_method="average"),
    cem_scm = get_att_cem(df, num_bins=5, est_method="scm"),
    cem_avg = get_att_cem(df, num_bins=5, est_method="average"),
    onenn = get_att_1nn(df, num_bins=5),
    
    # tmle1 = get_att_tmle(df, covs=c(X1,X2),
    #                      Q.SL.library = SL.library1,
    #                      g.SL.library = SL.library1),
    # tmle2 = get_att_tmle(df, covs=c(X1,X2),
    #                      Q.SL.library = SL.library3Q,
    #                      g.SL.library = SL.library3g),
    # tmle3 = get_att_tmle(df, covs=c(X1,X2,X3),
    #                      Q.SL.library = SL.library3Q,
    #                      g.SL.library = SL.library3g),

    # aipw1 = get_att_aipw(df, covs=c(X1,X2,X3),
    #                      Q.SL.library = SL.library1,
    #                      g.SL.library = SL.library1),
    # aipw2 = get_att_aipw(df, covs=c(X1,X2),
    #                      Q.SL.library = SL.library2,
    #                      g.SL.library = SL.library2),
    # aipw3 = get_att_aipw(df, covs=c(X1,X2),
    #                      Q.SL.library = SL.library2,
    #                      g.SL.library = SL.library2)
  )
  
  res %>% 
    pivot_longer(-c(runid:true_ATT))
  
  FNAME <- "sim_adversarial_results/toy_constant-tx-effect.csv"
  if (file.exists(FNAME)) {
    write_csv(res, FNAME, append=T)
  } else {
    write_csv(res, FNAME)
  }
}
toc()



# NOTES:
#  - test_full2: tx_effect = function(X1, X2) {(X1-0.5)^2+(X2-0.5)^2}
#  - test_full3: tx_effect = function(X1, X2) {X1+X2}
#  - test_full4: tx_effect = function(X1, X2) {(X1+X2)*2}, less sd
#  - test_ps: check that ps works


# analyze results ---------------------------------------------------------

# plot sample data
set.seed(1)
# dat <- gen_df_advhet(nc=1000, nt=50)
dat <- gen_df_full(nc=1000, nt=100)
vlines <- seq(min(dat$X1),max(dat$X1),length.out=6)
hlines <- seq(min(dat$X2),max(dat$X2),length.out=6)
dat %>% 
  ggplot(aes(x=X1, y=X2, color=Y)) +
  geom_point(aes(pch=Z)) +
  scale_color_continuous(low="blue", high="orange") +
  geom_vline(xintercept=vlines, lty="dotted") +
  geom_hline(yintercept=hlines, lty="dotted")




total_res <- read_csv("sim_adversarial_results/test_full4.csv")
# show results

total_res %>% 
  pivot_longer(-runid) %>% 
  group_by(name) %>% 
  summarize(rmse = sqrt(mean((value-tx_effect)^2)),
            bias = mean(value-tx_effect),
            sd   = sd(value)) %>% 
  rename(method = name) %>% 
  # pivot_longer(c(rmse, bias, sd)) %>% 
  ggplot(aes(x=method)) +
  geom_col(aes(y=rmse)) +
  coord_flip()


total_res %>% 
  pivot_longer(-(runid:true_ATT)) %>%
  filter(!name %in% c("ps_lm", "ps_bart")) %>% 
  group_by(nc,nt,eps_sd,name) %>%
  summarize(rmse = sqrt(mean((value-true_ATT)^2)),
            bias = mean(value-true_ATT),
            sd   = sd(value)) %>% 
  rename(method = name) %>% 
  pivot_longer(c(rmse, bias, sd)) %>%
  ggplot(aes(x=method)) +
  geom_col(aes(y=value)) +
  facet_wrap(~name) +
  coord_flip()

# CEM does worse now, since it suffers some bias (due to dropping units with extreme tx effects)
# but 1nn still completely nails it! As does BART, but that's okay.

