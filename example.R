
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
# df <- df

# use experimental treateds and nonexperimental controls
df <- df %>%
  filter((dataset == "NSWD" & Z==T) | (dataset=="CPS" & Z==F))


# EDA ---------------------------------------------------------------------

df %>%
  group_by(dataset, Z) %>%
  summarize(across(starts_with("X"),
                   mean))

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
  theme(legend.position = "bottom")
ggsave("writeup/figures/lalonde_dist.png", height=3, width=4)

# compare ESS
ess_plot(feasible)
ggsave("writeup/figures/lalonde_ess.png", height=3, width=4)



# (optional) visualize density --------------------------------------------

require(latex2exp)
dm <- data.frame(t(as.matrix(attr(calada_scm, "dm_uncapped"))))
# rank each column
dm_col_sorted <- apply(dist_matrix, 2, sort)
dist_to_plot <- dm_col_sorted[1:3,] # choose top 3
tibble(d = as.numeric(as.matrix(dist_to_plot))) %>%
  filter(d < 1000) %>%
  ggplot(aes(d)) +
  geom_histogram(color="black", binwidth=0.3) +
  theme_classic() +
  labs(y = "Count",
       x = TeX("$d(X_t, X_j)$"))
ggsave("writeup/figures/lalonde_calselect_top_3.png", height=2.5, width=5)


# tibble(d = as.numeric(attr(calada_scm, "dm_uncapped"))) %>%
#   filter(d < 1000) %>%
#   ggplot(aes(d)) +
#   geom_histogram(color="black", binwidth=0.5) +
#   theme_classic() +
#   labs(y = "Count",
#        x = TeX("$d(X_t, X_j)$"))
# ggsave("writeup/figures/lalonde_calselect.png", height=2.5, width=5)



# check some units --------------------------------------------------------

calada_scm %>%
  left_join(attr(calada_scm, "adacalipers") %>%
              rename(subclass = id),
            by = "subclass") %>%
  arrange(desc(adacal)) %>%
  filter(adacal > 1, Z==1) %>%
  select(id, X5:X8, adacal)


# estimate-estimand tradeoff diagnostics ----------------------------------

# maximum caliper vs. # co units added
set.seed(1)
# require(patchwork)
# satt_plot(calada_scm, B=NA) +
#   geom_hline(yintercept=0, lty="dotted")
# satt_plot2(calada_scm, B=100) +
#   geom_hline(yintercept=0, lty="dotted")
# foo <- satt_plot4(calada_scm, B=1000)
foo <- satt_plot4(calada_scm, B=NA)
foo +
  geom_hline(yintercept=0, lty="dotted") +
  theme(legend.direction="horizontal",
        legend.position = c(0.5, 0.85),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  labs(y = "Cumulative ATT Estimate",
       x = "Total number of treated units used",
       color = "Maximum \ncaliper \nsize used    ")+
  ylim(c(0,2500))
ggsave("writeup/figures/lalonde_att.png", height=3, width=6)
# [1] "FSATT: (339.098, 2851.353)"
# [1] "SATT: (76.115, 2612.28)"


# Love plot vs. # co units added
set.seed(2)
love_plot(calada_scm, covs=paste0("X",5:8), B=NA) +
  scale_color_manual(values = wesanderson::wes_palette("Zissou1", 5)[c(5,3,2,1)]) +
  guides(color=F) +
  scale_x_continuous(breaks = c(174, 177, 180, 183))
ggsave("writeup/figures/lalonde_love.png", height=3, width=4)


love_labs <- c("Age", "Years of Education", "Earnings (1974)", "Earnings (1975)")
names(love_labs) <- paste0("X", 5:8)
tmp <- attr(calada_scm,  "adacalipers")
love_plot2(calada_scm, covs = paste0("X", 5:8)) +
  theme(legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  scale_color_manual(values = wesanderson::wes_palette("Zissou1", 5)[c(1,5)]) +
  facet_wrap(~name,
             scales="free",
             labeller = labeller(name = love_labs))
ggsave("writeup/figures/lalonde_love2.png", height=3, width=6)



### naive bootstrap

boot_naive_res <- boot_naive(df,
                             caliper = CALIPER,
                             metric = METRIC,
                             cal_method = CAL_METHOD,
                             dist_scaling = DIST_SCALING,
                             B = 50)

write_csv(boot_naive_res, file="sim_canonical_results/boot_naive_res.csv")
