
source("./scripts/datagen/gen_six_points.R")
source("./R/wrappers.R")
source("./R/estimate.R")

### Step 1: generate data
df_six_points <-
  gen_six_points(tau=function(x1,x2) 0,
                 scale_Y0 = 10)
# we tune hardly to see how bad balancing


# Plot
ggplot(df_six_points,
       aes(x=X1,
           y=X2,
           pch=as.factor(Z))) +
  geom_point()

# True ATT = Y1 - Y0
true_att <- df_six_points %>%
  filter(Z==1) %>%
  mutate(eff= Y1-Y0) %>%
  group_by(Z) %>%
  summarise(mean_eff = mean(eff)) %>%
  pull(mean_eff)
true_att

# estimated ATT if you match to local controls, i.e.
#   what CSM is supposed to do
good_wt <- c(1,1,0,0,1,1)
(good_est_att <-
    get_est_att_from_wt(df=df_six_points,
                        input_wt=good_wt))

# estimated ATT if you match to all controls
# i.e. what
# weighted sum of
# wt = 1 for 1st, 2nd units (trts)
# wt = 2/3 for 3rd 4th unit (bad controls)
# wt = 1/3 for 5th, 6th unit (good controls)
# weighted sum of Y in 1st, 2nd - weighted sum of Y in 3rd-6th units
bad_wt <- c(1,1,2/3,2/3,1/3,1/3)
(bad_est_att <- get_est_att_from_wt(df=df_six_points,
                    input_wt=bad_wt))

##
covs <- c("X1", "X2")
zform1 <- as.formula(paste0("Z ~ ", paste0(covs, collapse="+")))
(bal_est_att = get_att_bal(df_six_points,
                   zform1,
                   rep(0.01, length(covs))))

expect_true(bad_est_att, bal_est_att, tolerance=0.0001)

# okay so bad_est_att and good_est_att does not
#   differ that much. This is really a bad news
# Reason this is bad:
# Solution: increase f_0(x) by scaling up



# Check that the center of
#  (c1 * 2, c2*2, c3, c4)
#  is still the center of the treated
bad_wt <- c(1,1,2/3,2/3,1/3,1/3)
good_wt <- c(1,1,0,0,1,1)
wt_test <-good_wt
wtd_centers_df <- df_six_points %>%
  cbind(wt=wt_test) %>%
  group_by(Z) %>%
  summarise(wtd_X1 = weighted.mean(X1, wt),
            wtd_X2 = weighted.mean(X2, wt))
wtd_centers_df
# Result (wtd_X1, wtd_X2) are close in trt and control group

### Step 2: make it work on CSM
# source("R/matching.R")
cal_matches <- get_cal_matches(df_six_points,
                metric = c("maximum"),
                caliper = 1,
                cal_method = "adaptive",
                est_method = "average",
                return = "sc_units",
                dist_scaling = 1)

### Step 2.5: Choose Mtd_Comp, and make it work
covs <- c("X1", "X2")
zform2 <- as.formula(paste0("Z ~ (", paste0(covs, collapse="+"), ")^2"))

# m_bal <- optweight(form=zform1,
#                    data = df_six_points,
#                    tols = rep(0.01, length(covs)),
#                    estimand = "ATT")
# m_bal$weights

bal2 = get_att_bal(df, zform2,
                   rep(0.1, length(covs) + choose(length(covs), 2)))


### Step 3: compute metrics of success
# Candidates:
# 0. matched units for t1
#   Expectation:
#     CSM: C3 only
#     Mtd_Comp: C3, C1, C2
# a. average over t1, t2, cov distance between trt and
#   Name: averaged matched controls.
#   Expectation: CSM < Mtd_Comp, because CSM does local matching
# b. average over t1, t2, difference between Y_t(0)- \hat Y_c(0)
#   Name:
#   Expectation:
