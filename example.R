
# example run

require(tidyverse)

source("R/distance.R")
source("R/sc.R")
source("R/matching.R")
source("R/estimate.R")
source("R/inference.R")
source("R/diagnostic_plots.R")


# NAME ASSUMPTIONS:
#  - covariates start with "X", and no other columns start with "X"
#  - outcomes start with "Y"
#  - treatment is "Z"

# DATA ASSUMPTIONS:
#  - exact match on logicals/categoricals (categoricals one-hot encoded)
#     - NOTE: all logicals/categoricals need to be T/F, not 0/1,
#       to distinguish that they should be exactly matched
#  - approximately match on numerics

# DISTANCE ASSUMPTIONS:
#  - default scaling is 1/sd(x)
#  - (default categorical scaling is 1000, i.e., some big number)


# load data ---------------------------------------------------------------

lalonde <- read_csv("lalonde.csv") %>% 
  rename(Z = treat, Y = re78) %>%
  mutate(id = 1:n(),
         dataset = c(rep("NSWD", 185+260),
                     rep("CPS", 15992)),
         across(c(black, hispanic, married, nodegree), as.logical)) %>% 
  select(id, Z, Y,
         black, hispanic, married, nodegree, 
         age, education, re74, re75,
         dataset)

# rename variables for my matching methods
df <- lalonde %>%
  rename(X1=black, X2=hispanic, 
         X3=married, X4=nodegree, 
         X5=age, X6=education,
         X7=re74, X8=re75)

# # for now, use only experimental data
# df <- df %>%
#   filter(dataset == "NSWD")

# # use experimental treateds and all controls
# df <- df %>% 
#   filter(Z==T | (dataset=="CPS" & Z==F))

# use experimental treateds and nonexperimental controls
df <- df %>% 
  filter((dataset == "NSWD" & Z==T) | (dataset=="CPS" & Z==F))


# EDA ---------------------------------------------------------------------

if (F) {
  # NSWD is...
  #  - younger, less married
  #  - less educated, fewer degrees
  #  - much more black
  #  - much lower earnings
  df %>% 
    group_by(dataset) %>% 
    summarize(across(X1:X8, mean))
  
  # NSWD has...
  #  - less variance in age/educ/income
  df %>% 
    group_by(dataset) %>% 
    summarize(across(X1:X8, sd))
  
  # visualize covariate distributions
  df %>% 
    group_by(dataset) %>% 
    pivot_longer(X1:X8) %>% 
    ggplot(aes(x=value, color=dataset)) +
    geom_density() +
    facet_wrap(~name, scales="free")
}


# set matching settings ---------------------------------------------------

CALIPER <- 1
METRIC <- "maximum"   # "maximum", "euclidean", "manhattan"
CAL_METHOD <- "ada"   # "ada", "cem"
DIST_SCALING <- tibble(
  X1 = 1000,
  X2 = 1000,
  X3 = 1000,
  X4 = 1000,
  X5 = 1/3,
  X6 = 1,
  X7 = 1/5000,
  X8 = 1/5000
)

# DIST_SCALING <- tibble(
#   X1 = 1/5,
#   X2 = 1/4,
#   X3 = 1000,
#   X4 = 1000,
#   X5 = 1000,
#   X6 = 1000,
#   X7 = 1/10000,
#   X8 = 1/10000
# )
# 
# DIST_SCALING <- tibble(
#   X1 = 1/10,
#   X2 = 1/4,
#   X3 = 1000,
#   X4 = 1000,
#   X5 = 1000,
#   X6 = 1000,
#   X7 = 1/20000,
#   X8 = 1/20000
# )


# run matching ------------------------------------------------------------

calada_scm <- df %>%
  get_cal_matches(caliper = CALIPER,
                  metric = METRIC,
                  cal_method = CAL_METHOD,
                  dist_scaling = DIST_SCALING,
                  est_method = "scm",
                  return = "sc_units",
                  knn = 10,         # for ada
                  num_bins = 5,     # for cem
                  wider = F)        # for cem
get_att_ests(calada_scm)

# get_cem_matches(df, num_bins=5, method="average", return="sc_units")



# FSATT results -----------------------------------------------------------

feasible <- attr(calada_scm, "scweights") %>% 
  bind_rows() %>% 
  filter(subclass %in% attr(calada_scm, "feasible_subclasses"))
get_att_ests(feasible)


# compare distances between units represented by average and sc weights
# scm_vs_avg_plot(feasible, DIST_SCALING, METRIC)
dist_density_plot(feasible, DIST_SCALING, METRIC) +
  theme(legend.position = c(0.82,0.68))
ggsave("writeup/figures/lalonde_dist.png", height=3, width=5)

# compare ESS
ess_plot(feasible)
ggsave("writeup/figures/lalonde_ess.png", height=3, width=5)


# estimate-estimand tradeoff diagnostics ----------------------------------

# maximum caliper vs. # co units added
require(patchwork)
set.seed(1)
satt_plot(calada_scm, B=100) +
  geom_hline(yintercept=0, lty="dotted")
ggsave("writeup/figures/lalonde_att.png", height=4, width=5)


# Love plot vs. # co units added
set.seed(2)
love_plot(calada_scm, covs=paste0("X",5:8), B=100) +
  scale_color_manual(values = wesanderson::wes_palette("Zissou1", 5)[c(5,3,2,1)]) +
  guides(color=F)
ggsave("writeup/figures/lalonde_love.png", height=3, width=5)






### naive bootstrap

boot_naive_res <- boot_naive(df, 
                             caliper = CALIPER,
                             metric = METRIC,
                             cal_method = CAL_METHOD,
                             dist_scaling = DIST_SCALING,
                             B = 50)


