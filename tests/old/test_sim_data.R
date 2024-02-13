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

test_matches_from_otsu <- function(){
  df_otsu <- generate_one_otsu()
  otsu_matched<- get_NN_matches(df_otsu)
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
    
  }


# library(MatchIt)
# m.out <- matchit(Z ~ X1 + X2, 
#                  data = df_otsu, 
#                  method = "nearest", 
#                  replace = TRUE,
#                  mahvars = ~ X1 + X2,
#                  ratio = 8)
# 
# # Obtain matched dataset
# df_matched <- match.data(m.out)
# 
# df_matched_ids <- data.frame(m.out$match.matrix)
# 
# id <- subclass <- weights <- c()
# 
# for(i in 1:nrow(df_matched_ids)) {
#   # Concatenate row name (treated id) with matched control ids
#   # i <- 1
#   row_ids <- c(as.integer(rownames(df_matched_ids)[i]), 
#                as.integer(unlist(df_matched_ids[i, ])))
#   id <- c(id, row_ids)
#   # Add subclass
#   M <- length(row_ids) - 1
#   subclass <- c(subclass, rep(i, M + 1))
#   weights <- c(weights, c(1, rep(1/M, M )) )
# }
# 
# # Create the new dataframe
# new_df <- data.frame(id = id, 
#                      subclass = subclass,
#                      weights)
# 
# # Display the first few rows of the new dataframe
# head(new_df)
# 
# preds_matching_otsu <- 
#   left_join(new_df, df_otsu, by = "id")
