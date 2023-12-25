setwd("~/Dropbox (Harvard University)/Xiang_Luke/scmatch2")
library(tidyverse)
require(mvtnorm)
source("./R/sim_data.R")

test_gen_df_otsu <- function(){
    df_dgp <- generate_one_toy()
    df_otsu <- generate_one_otsu()
    ggplot(df_otsu, aes(x=X1,y=X2))+
      geom_point()
    
}
source("./R/sim_data.R")
test_get_df_scaling_from_otsu <- function(){
  dgp_obj <- 
    get_df_scaling_from_dgp_name(dgp_name="otsu",
                                 kang_true=F)
  list2env(dgp_obj,envir = environment())
}

test_matches_and_debiased_residuals_from_otsu <-
  function(){
    # First get one from toy. See the return format
    dgp_obj <- 
      get_df_scaling_from_dgp_name(dgp_name="toy",
                                   kang_true=F)
    list2env(dgp_obj,envir = environment())
    match_obj_toy <-
      get_matches_and_debiased_residuals(
      "toy", df_dgp, 
      dist_scaling, "linear",n_split=2)
    list2env(match_obj_toy, 
             envir = environment())
    preds_csm
    
    tilde_tau
    tilde_tau_resids
    # Then get from
    
  }

df_otsu <- generate_one_otsu()
head(df_otsu)
library(MatchIt)
m.out <- matchit(Z ~ X1 + X2, 
                 data = df_otsu, 
                 method = "nearest", 
                 replace = TRUE,
                 mahvars = ~ X1 + X2,
                 ratio = 8)

# Obtain matched dataset
df_matched <- match.data(m.out)

df_matched_ids <- data.frame(m.out$match.matrix)



