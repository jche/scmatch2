
library( tidyverse )
library( CSM )
source( here::here( "scripts/datagen/gen_uniform.R" ) )

xlim <- ylim<-  c(0, 2.25)
uniform_df <- gen_uniform(n = 500,
                          p_trt = 0.5,
                          xlim = xlim,
                          ylim = ylim)

ggplot(uniform_df,
       aes(x=X1,
           y=X2,
           color=as.factor(Z),
           pch=as.factor(Z))) +
  geom_point()

### methods:
## - CSM: average weight
## - CEM: radius
DIST_SCALING = 1 # Keep this fixed


x_range <- (xlim[2] - xlim[1])



res <- NULL
R = 100
for (num_bins in c(4, 1)){
  num_bins = 1
  CEM_size = x_range / num_bins
  caliper = CEM_size / 2

  CSM_ests <- CEM_ests <- numeric(R)
  set.seed(119+num_bins)
  for (r in 1:R){
    print(r)
    uniform_df <- gen_uniform(n = 500,
                              p_trt = 0.5,
                              xlim = xlim,
                              ylim = ylim)
    cal_matches <-
      get_cal_matches(uniform_df,
                      covs=c("X1","X2"),
                      treatment = "Z",
                      metric = "maximum",
                      caliper = caliper,  # impt: caliper
                      cal_method = "adaptive", # not important
                      est_method = "scm", # impt: weighting method
                      return = "all",
                      dist_scaling = DIST_SCALING
      )
    cal_matches$Y <- cal_matches$Y0
    CSM_ests[r] <- CSM_est <- get_att_ests(cal_matches)


    ## CEM
    tmp <- MatchIt::matchit(as.formula(Z ~ X1 + X2),
                            data = uniform_df,
                            method = "cem",
                            estimand = "ATT",
                            s.weights = NULL,
                            verbose = FALSE,
                            grouping = NULL,   # exact match factor covariates
                            cutpoints = num_bins,
                            k2k = F)
    m.data3 <- MatchIt::match.data(tmp)
    m.data3$Y <- m.data3$Y0
    CEM_ests[r] <- CEM_est <-
      get_att_ests(m.data3)

  }
  res_curr_bin_size <- data.frame(
    run = 1:R,
    CSM_ests=CSM_ests,
    CEM_ests=CEM_ests,
    CEM_num_bins = num_bins,
    CSM_radius = caliper
  )

  res<- rbind(res, res_curr_bin_size)
}




saveRDS(res,
        file="data/outputs/CSM-vs-CEM-results/CSM-vs-CEM-uniform.rds")

res <- readRDS(file="data/outputs/CSM-vs-CEM-results/CSM-vs-CEM-uniform.rds")

res_summary <- res %>%
  group_by(CEM_num_bins, CSM_radius) %>%
  summarise(
    n_runs = n(),
    sd_CSM = sd(CSM_ests),
         sd_CEM = sd(CEM_ests),
         abs_bias_CSM = mean(abs(CSM_ests)),
         abs_bias_CEM = mean(abs(CEM_ests))
  )


## A. What the above said to inference is:
#   With this Y(0) generation, CSM with caliper
#   size 0.28125, the SD(bias) is 0.0002766.
# Note the largest Y(0) at c(1.5,1.5) is 1.8, so
#   the SD(bias) where we randomized treatment
#   and control covariates are small.
# The question: does whole uniform sampling
#   reduces the SD(bias)?
# The answer here: yes. There is also the component
#   of ratio
#   The Y(0) generation:
MU = c(1.5,1.5)
SIG = matrix(c(1,0.5,0.5,1), nrow=2)
x1_test = 1.5; x2_test = 1.5
mvtnorm::dmvnorm(x = c(x1_test, x2_test),
                 mean = MU,
                 sigma = SIG)


# In existing inference result, Y(0) in the toy
#   example is made by the following function
20 * mvtnorm::dmvnorm(x = c(0.5, 0.5),
                 mean = c(0.5,0.5),
                 sigma = matrix(c(1,0.8,0.8,1), nrow=2))
# The radius is 6/1 = 6, so it is not matching
# Is it reasonable that we have expected

## B. What the above said to CSM vs CEM is:
#  Compared to CEM, CSM has slightly larger bias
#     but smaller variance.



# cem_scm = get_att_cem(uniform_df,
#                       num_bins=num_bins,
#                       est_method="scm",
#                       estimand = "CEM-ATT")
# cem_avg = get_att_cem(uniform_df,
#                       num_bins=num_bins,
#                       est_method="average",
#                       estimand = "CEM-ATT")
#
#
# Play: vary scaling factor, i.e., radius.
