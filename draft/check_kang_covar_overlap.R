setwd("~/Dropbox (Harvard University)/Xiang_Luke/scmatch2")
library(tidyverse)
require(mvtnorm)
source("R/distance.R")
source("R/sc.R")
source("R/matching.R")
source("R/estimate.R")
source("R/inference.R")
source("R/sim_data.R")
source("R/wrappers.R")
source("R/utils.R")
source("R/bootstrap.R")

df_kang <- gen_df_kang(n = 1000)

ggplot(data=df_kang, aes(x=V1, y=V2)) +
  geom_point(aes(color=Z))

