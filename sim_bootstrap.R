
# check validity of bootstrap

require(tidyverse)
require(mvtnorm)

require(tictoc)

source("R/distance.R")
source("R/sc.R")
source("R/matching.R")
source("R/estimate.R")
source("R/inference.R")
source("R/sim_data.R")
source("R/wrappers.R")
source("R/utils.R")

# parallelize
if (T) {
  library(foreach)
  library(doParallel)
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
}



# function to get bootstrap results ---------------------------------------

# feels like original approach was correct, according to second line of (3.1)
# new approach uses first line of (3.1)
#  --> idea is just that you can get att from sc units, OR from weighted original units

res_fun <- function(d, dist_scaling, name, att0, B=250) {

  # # get CSM results
  # preds_csm <- get_cal_matches(
  #   df = d,
  #   metric = "maximum",
  #   dist_scaling = dist_scaling,
  #   cal_method = "fixed",
  #   est_method = "scm",
  #   return = "all",   # maybe just "agg_co_units"...?
  #   knn = 25)
  #
  # # bootstrap FSATT
  # boot_fsatt <- attr(preds_csm, "scweights") %>%
  #   bind_rows() %>%
  #   filter(subclass %in% attr(preds_csm, "feasible_subclasses")) %>%
  #   agg_co_units() %>%
  #   boot_bayesian(B=B)

  # get CSM results
  preds_csm <- get_cal_matches(
    df = d,
    metric = "maximum",
    dist_scaling = dist_scaling,
    cal_method = "fixed",
    est_method = "scm",
    return = "sc_units",
    knn = 25)   # NOTE: knn does nothing when cal_method = "fixed"...

  if (F) {
    # check distribution of aggregate co weights
    foo <- get_cal_matches(
      df = d,
      metric = "maximum",
      dist_scaling = dist_scaling,
      cal_method = "fixed",
      est_method = "scm",
      return = "agg_co_units",
      knn = 25)
    foo %>%
      filter(!Z) %>%
      pull(weights) %>%
      hist()

    # check balance
    preds_csm %>%
      pivot_longer(starts_with("X")) %>%
      ggplot(aes(x=value, color=Z)) +
      geom_density() +
      facet_wrap(~name, scales="free")
  }

  boot_fsatt <- preds_csm %>%
    boot_bayesian2(B=B)

  if (att0) {
    att <- 0
  } else {
    # grab true att among feasible units
    att <- d %>%
      # filter((!Z) | id %in% attr(preds_csm, "feasible_units")) %>%
      filter(Z & !(id %in% attr(preds_csm, "unmatched_units"))) %>%
      summarize(att = mean(Y1-Y0)) %>%
      pull(att)
  }

  tibble(
    id = i,
    data = name,
    att = att,
    att_hat = get_att_ests(preds_csm),
    sd = sd(boot_fsatt),
    q025 = quantile(boot_fsatt, 0.025),
    q975 = quantile(boot_fsatt, 0.975)
  ) %>%
    mutate(covered1 = att_hat-1.96*sd < att & att < att_hat+1.96*sd,
           covered2 = q025 < att & att < q975,
           length1 = 1.96*2*sd,
           length2 = q975-q025)
}

res_fun_wild <- function(d, dist_scaling, name, att0, B=250) {

  # get CSM results
  preds_csm <- get_cal_matches(
    df = d,
    metric = "maximum",
    dist_scaling = dist_scaling,
    cal_method = "fixed",
    est_method = "scm",
    return = "sc_units",
    knn = 25)   # NOTE: knn does nothing when cal_method = "fixed"...

  boot_fsatt <- preds_csm %>%
    boot_bayesian2_wild(B=B)

  if (att0) {
    att <- 0
  } else {
    # grab true att among feasible units
    att <- d %>%
      # filter((!Z) | id %in% attr(preds_csm, "feasible_units")) %>%
      filter(Z & !(id %in% attr(preds_csm, "unmatched_units"))) %>%
      summarize(att = mean(Y1-Y0)) %>%
      pull(att)
  }

  tibble(
    id = i,
    data = name,
    att = att,
    att_hat = get_att_ests(preds_csm),
    sd = sd(boot_fsatt)
  ) %>%
    mutate(covered = att_hat-1.96*sd < att & att < att_hat+1.96*sd,
           length = 1.96*2*sd)
}


# super wide...
res_fun_naive <- function(d, dist_scaling, name, att0, B=50) {
  boot_samps <- boot_naive(
    d,
    caliper = 1,
    metric = "maximum",
    dist_scaling = dist_scaling,
    cal_method = "fixed",
    est_method = "scm",
    knn = 25,
    B = B)

  browser()
}


res_fun_redux <- function(d, dist_scaling, name, att0, B=250) {

  # get CSM results, aggregated co units
  preds_csm <- get_cal_matches(
    df = d,
    metric = "maximum",
    dist_scaling = dist_scaling,
    cal_method = "fixed",
    est_method = "scm",
    return = "agg_co_units",
    knn = 25)   # NOTE: knn does nothing when cal_method = "fixed"...

  boot_fsatt <- preds_csm %>%
    boot_bayesian_finalattempt(B=B)

  if (F) {
    hist(boot_fsatt)
    get_att_ests(preds_csm)

    preds_csm %>%
      filter(Z) %>%
      summarize(att = mean(Y1-Y0))
  }

  if (att0) {
    att <- 0
  } else {
    # grab true att among feasible units
    att <- d %>%
      # filter((!Z) | id %in% attr(preds_csm, "feasible_units")) %>%
      filter(Z & !(id %in% attr(preds_csm, "unmatched_units"))) %>%
      summarize(att = mean(Y1-Y0)) %>%
      pull(att)
  }

  tibble(
    id = i,
    data = name,
    att = att,
    att_hat = get_att_ests(preds_csm),
    q025 = quantile(boot_fsatt, 0.025),
    q975 = quantile(boot_fsatt, 0.975)
  ) %>%
    mutate(covered = att_hat-q975 < att & att < att_hat-q025)
    # mutate(covered = q025 < att & att < q975)
    # mutate(covered = att_hat-q975 < att & att < att_hat-q025)
}



# run simulations ---------------------------------------------------------

# repeatedly run for all combinations of pars
res <- foreach(
  i=1:100,
  .packages = c("tidyverse", "mvtnorm"),
  .combine=rbind) %dopar% {
    source("R/distance.R")
    source("R/sc.R")
    source("R/matching.R")
    source("R/estimate.R")
    source("R/inference.R")
    source("R/sim_data.R")
    source("R/wrappers.R")
    source("R/utils.R")
    
    # for (i in 1:250) {
    
    # toy example
    df_toy <- gen_df_adv(
      nc=500, 
      nt=100, 
      f0_sd = 0.5,
      tx_effect_fun = function(X1, X2) {3*X1+3*X2},
      f0_fun = function(x,y) {
        matrix(c(x,y), ncol=2) %>%
          dmvnorm(mean = c(0.5,0.5),
                  sigma = matrix(c(1,0.8,0.8,1), nrow=2)) * 20   # multiply for more slope!
      })
    dist_scaling_toy <- df_toy %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) 6 / (max(x) - min(x))
                         else 1000
                       }))
    
    # kang
    # set.seed(2)   # for example of bad results
    df_kang <- gen_df_kang(n = 1000)
    dist_scaling_kang <- df_kang %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) 5 / (max(x) - min(x))
                         else 1000
                       }))
    
    # hain
    df_hain <- gen_df_hain(
      nc=250, 
      nt=50, 
      sigma_e = "n100",   # high overlap condition
      outcome = "nl2",
      sigma_y = 1,
      ATE = 0)
    dist_scaling_hain <- df_hain %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) 5 / (max(x) - min(x))
                         else 1000
                       }))
    
    # acic
    df_acic <- gen_df_acic(
      model.trt="step", 
      root.trt=0.35, 
      overlap.trt="full",
      model.rsp="step", 
      alignment=0.75, 
      te.hetero="high",
      random.seed=i,        # NOTE: set random seed!
      n=1000, 
      p=10
    )
    dist_scaling_acic <- df_acic %>%
      summarize(across(starts_with("X"),
                       function(x) {
                         if (is.numeric(x)) 5 / (max(x) - min(x))
                         else 1000
                       }))
    
    
    # get results
    rf <- res_fun_redux
    res_toy <- rf(df_toy, dist_scaling_toy, name="toy", att0=F)
    res_kang <- rf(df_kang, dist_scaling_kang, name="kang", att0=T)
    res_hain <- rf(df_hain, dist_scaling_hain, name="hain", att0=T)
    res_acic <- rf(df_acic, dist_scaling_acic, name="acic", att0=F)
    res <- bind_rows(res_toy, res_kang, res_hain, res_acic)
    
    FNAME <- "bootstrap_results/bootstrap_results_finalattempt.csv"
    if (file.exists(FNAME)) {
      write_csv(res, FNAME, append=T)
    } else {
      write_csv(res, FNAME)
    }
    
    res
  }

stopCluster(cl)



# analyze results ---------------------------------------------------------

# for "fixed" intervals, we were computing acic att wrong
# for "fixed2" intervals, acic/kang coverage still very bad...
#  - Q: why are kang point estimates sometimes so bad,
#       even with so many control matches per tx unit?

res <- read_csv("bootstrap_results/bootstrap_results_fixed2.csv")

res %>% 
  group_by(data) %>% 
  summarize(coverage = mean(covered),
            avg_sd = mean(sd),
            avg_len = mean(length),
            avg_bias = mean(abs(att-att_hat)))

# check bias of estimates
res %>% 
  group_by(data) %>% 
  pivot_longer(c(att, att_hat)) %>% 
  ggplot() +
  geom_density(aes(x=value, color=name)) +
  facet_wrap(~data, scales="free")

# check distribution of biases
res %>% 
  ggplot() +
  geom_density(aes(x=att_hat-att)) +
  facet_wrap(~data, scales="free")

# check how long CIs should be
res %>% 
  group_by(data) %>% 
  summarize(q5 = quantile(att_hat-att, 0.025),
            q95 = quantile(att_hat-att, 0.975)) %>% 
  mutate(len = q95-q5)




res %>% 
  group_by(data) %>% 
  summarize(sd_att = sd(att),
            sd_hat = sd(att_hat))









res <- read_csv("bootstrap_results/bootstrap_results_q.csv")
res <- read_csv("bootstrap_results/bootstrap_results_finalfornow.csv")


res %>% 
  group_by(data) %>% 
  summarize(coverage1 = mean(covered1),
            coverage2 = mean(covered2),
            avg_len1 = mean(length1),
            avg_len2 = mean(length2))



res <- read_csv("bootstrap_results/bootstrap_results_redux.csv")

res %>% 
  group_by(data) %>% 
  summarize(coverage = mean(covered))

res %>% 
  group_by(data) %>% 
  summarize(avg_len = mean(q975-q025))




res <- read_csv("bootstrap_results/bootstrap_results_redux2.csv")

res %>% 
  group_by(data) %>% 
  summarize(coverage = mean(covered))

res %>% 
  group_by(data) %>% 
  summarize(avg_len = mean(q975-q025))





res <- read_csv("bootstrap_results/bootstrap_results_redux3.csv")









