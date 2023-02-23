
# example run

require(tidyverse)

source("R/distance.R")
source("R/sc.R")
source("R/matching.R")
source("R/estimate.R")


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
         age, education, black, hispanic, married, nodegree, re74, re75,
         dataset)

# rename variables for my matching methods
df <- lalonde %>%
  rename(X1=age, X2=education, 
         X3=black, X4=hispanic, 
         X5=married, X6=nodegree, 
         X7=re74, X8=re75)

# for now, use only experimental data
df <- df %>% 
  filter(dataset == "NSWD")


# run matching ------------------------------------------------------------

CALIPER <- 1
METRIC <- "manhattan"   # "maximum", "euclidean", "manhattan"
CAL_METHOD <- "ada"   # "ada", "cem"
DIST_SCALING <- tibble(
  X1 = 1/3,
  X2 = 1/2,
  X3 = 1000,
  X4 = 1000,
  X5 = 1000,
  X6 = 1000,
  X7 = 1/5000,
  X8 = 1/5000
)

calada_scm <- df %>%
  get_cal_matches(caliper = CALIPER,
                  metric = METRIC,
                  cal_method = CAL_METHOD,
                  dist_scaling = DIST_SCALING,
                  est_method = "scm",
                  return = "sc_units",
                  knn = 5,          # for ada
                  num_bins = 5,     # for cem
                  wider = F)        # for cem

attr(calada_scm, "scaling")
attr(calada_scm, "dm")[1:5,1:5]
attr(calada_scm, "adacalipers")
attr(calada_scm, "unmatched_units")

get_att_ests(calada_scm)


get_cem_matches(df, num_bins=5, method="average", return="sc_units")



# FSATT results -----------------------------------------------------------

# note: this section uses return="sc_units" output

# get feasible tx units
feasible_units <- attr(calada_scm, "adacalipers") %>% 
  filter(adacal <= CALIPER) %>% 
  pull(id)
feasible_subclasses <- calada_scm %>% 
  filter(id %in% feasible_units) %>% 
  pull(subclass)
calada_scm_feasible <- calada_scm %>% 
  filter(subclass %in% feasible_subclasses)

# check number of feasible units
length(feasible_units)

# check att estimate
get_att_ests(calada_scm_feasible)

# check distances bw each tx/sc pair
sc_dists <- gen_dm(calada_scm_feasible,
                   scaling = DIST_SCALING,
                   method = METRIC) %>% 
  diag()
calada_scm_feasible %>% 
  filter(Z==T) %>% 
  mutate(dist = sc_dists) %>% 
  ggplot(aes(x=dist)) +
  geom_density(linewidth=1) +
  geom_vline(xintercept=CALIPER, lty="dashed") +
  theme_classic() +
  labs(y = "",
       x = "Scaled distance between tx and sc unit") +
  theme_classic()

# compare distances to calipers
#  - see how well scm brings points down from x=y line
calada_scm_feasible %>% 
  filter(Z==T) %>% 
  mutate(dist = sc_dists) %>% 
  left_join(attr(calada_scm, "adacalipers"),
            by = "id") %>% 
  ggplot(aes(x=adacal, y=dist)) +
  geom_point() +
  geom_abline(lty="dotted") +
  theme_classic() +
  labs(x = "Adaptive caliper",
       y = "Dist. between tx and sc unit")


# TODO: more analysis on how feasible/infeasible subsamples differ

# "love plot" for scaled (feasible - infeasible) diffs
#  - uses scaled cov dists
scaled_diffs <- calada_scm %>% 
  filter(Z==T) %>% 
  mutate(feasible = subclass %in% feasible_subclasses) %>% 
  group_by(feasible) %>% 
  summarize(across(starts_with("X"), mean)) %>% 
  summarize(across(starts_with("X"), ~.x[2]-.x[1])) * 
  (DIST_SCALING %>% 
     mutate(across(everything(),
                   ~ifelse(.x==1000, 1, .x))))
scaled_diffs %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(x=name, y=value)) +
  geom_col() +
  geom_hline(yintercept=0, lty="dotted") +
  coord_flip() +
  theme_classic() +
  labs(y = "Scaled difference between feasible and infeasible subsamples",
       x = "Covariate")



# estimate-estimand tradeoff diagnostics ----------------------------------

# maximum caliper vs. # co units added
calada_scm %>% 
  filter(!subclass %in% feasible_subclasses) %>% 
  filter(Z==1) %>% 
  left_join(attr(calada_scm, "adacalipers"), by="id") %>% 
  arrange(adacal) %>% 
  mutate(order = 1:n()) %>% 
  ggplot(aes(x=order, y=adacal)) +
  geom_col() +
  theme_classic() +
  labs(y = "Adaptive caliper",
       x = "Unit #")

# ATT estimate vs. # co units added
calada_scm %>%
  left_join(attr(calada_scm, "adacalipers"), by="id") %>% 
  group_by(subclass) %>%
  summarize(adacal = last(adacal),
            tx = Y[2] - Y[1]) %>%
  arrange(adacal) %>%
  mutate(order = 1:n(),
         cum_avg = cumsum(tx) / order) %>% 
  slice(length(feasible_units):n()) %>% 
  ggplot(aes(x=order, y=cum_avg)) +
  geom_point() +
  geom_step(linewidth=1) +
  theme_classic() +
  labs(y = "Cumulative ATT Estimate",
       x = "Unit #")


# Love plot vs. # co units added
calada_scm %>%
  left_join(attr(calada_scm, "adacalipers"), by="id") %>% 
  group_by(subclass) %>%
  summarize(adacal = last(adacal),
            X1 = X1[2] - X1[1],
            X2 = X2[2] - X2[1],
            X7 = X7[2] - X7[1],
            X8 = X8[2] - X8[1]) %>%
  arrange(adacal) %>%
  mutate(order = 1:n(),
         Age = cumsum(X1) / order,
         Educ = cumsum(X2) / order,
         Income74 = cumsum(X7) / order,
         Income75 = cumsum(X8) / order) %>% 
  slice(length(feasible_units):n()) %>%
  pivot_longer(Age:Income75) %>% 
  ggplot(aes(x=order, y=value)) +
  geom_point(data=. %>% slice(1:4), 
             aes(color=name), size=2) +
  geom_step(aes(group=name, color=name), size=1.1) +
  expand_limits(y = 0) +
  geom_hline(yintercept=0, lty="dotted") +
  facet_wrap(~name, scales="free_y") +
  labs(y = "\n Covariate balance (tx-co)",
       x = "Maximum caliper size used",
       color = "Covariate") +
  scale_color_manual(values = wesanderson::wes_palette("Zissou1", 5)[c(1,2,3,5)]) +
  theme_classic()








