setwd("~/Dropbox (Harvard University)/Xiang_Luke/scmatch2")
library(tidyverse)
require(mvtnorm)
source("R/distance.R")
source("R/sc.R")
source("R/matching.R")
source("R/estimate.R")
source("R/inference.R")
source("R/sim_data.R")
source("scripts/wrappers.R")
source("R/utils.R")
source("R/bootstrap.R")

df_otsu <- generate_one_otsu()
library(MatchIt)
mNN <- matchit(Z ~ X1 + X2, 
               data = df_otsu,
               method = "nearest", 
               # replace = TRUE,
               mahvars = ~ X1 + X2,
               ratio = 8)
mNN
md <- match.data(mNN)

pair_ids <- levels(md$subclass)

#Unit IDs, split by pair membership
split_inds <- split(seq_len(nrow(md)), md$subclass)

cluster_boot_fun <- function(pairs, i) {
  
  n_ids <- length(pair_ids)
  # Sample with replacement on pair_ids, aka, treated units 
  i <- sample(1:n_ids,replace=T)
  # Since each treated is matched to M controls, 
  #   sampling in this way is sampling the whole dataset. 
  
  # Getting the resampled data. There might be repeated pairs
  ids <- unlist(split_inds[pair_ids[i]])
  
  #Subset md with block bootstrapped indices
  boot_md <- md[ids,]
  
  #Fit outcome model
  fit <- glm(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + 
                          X6 + X7 + X8 + X9),
             data = boot_md, weights = weights,
             family = quasibinomial())
  
  ## G-computation ##
  #Subset to treated units for ATT; skip for ATE
  md1 <- subset(boot_md, A == 1)
  
  #Estimated potential outcomes under treatment
  p1 <- predict(fit, type = "response",
                newdata = transform(md1, A = 1))
  Ep1 <- weighted.mean(p1, md1$weights)
  
  #Estimated potential outcomes under control
  p0 <- predict(fit, type = "response",
                newdata = transform(md1, A = 0))
  Ep0 <- weighted.mean(p0, md1$weights)
  
  #Risk ratio
  return(Ep1 / Ep0)
}
