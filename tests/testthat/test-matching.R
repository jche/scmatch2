
# test_df <-
#   data.frame(Z=c(1,0,0,0),
#              X=c(1,2,3,4))
#
source("./scripts/datagen/gen_six_points.R")
df_six_points <-
  gen_six_points(tau=function(x1,x2) 0,
                 scale_Y0 = 10)

res <- gen_matches(df=df_six_points,
                  covs = starts_with("X"),
                  treatment = "Z")


test_that("get_matched_co_from_dm should ")
