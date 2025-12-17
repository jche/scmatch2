# scripts/demo/example-lalonde.R
# example run

require(tidyverse)
library( CSM )


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
#  - (default categorical scaling is 1000, i.e., some big number, which forces exact matching on categoricals)


# load data ---------------------------------------------------------------
LOAD_OPTION <- "lalonde_w_cps_cos"
# option 1: "lalonde_only": use only experimental data
# option 2: "lalonde_w_cps": use experimental treats and all controls
# option 3: "lalonde_w_cps_cos": use experimental treateds and nonexperimental controls
load(file = paste0("data/inputs/", LOAD_OPTION,".RData"))


ggplot(lalonde_df,
       aes(x = married, y=re74 ))+
  geom_boxplot()+
  facet_wrap(.~Z)

# set matching settings ---------------------------------------------------

CALIPER <- 1
METRIC <- "maximum"   # "maximum", "euclidean", "manhattan"
CAL_METHOD <- "ada"   # "ada", "cem"

lalonde_df_renamed <- lalonde_df %>%
  rename(X1=black, X2=hispanic,
         X3=married, X4=nodegree,
         X5=age, X6=education,
         X7=re74, X8=re75)

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
#   black = 1000,
#   hispanic = 1000,
#   married = 1000,
#   nodegree = 1000,
#   age = 1/3,
#   education = 1,
#   re74 = 1/5000,
#   re75 = 1/5000
# )


# Number of treated units
sum( lalonde_df_renamed$Z )
sum( lalonde_df_renamed$Z ) * 2

# run matching
# TODO/NOTE: LWM reduced caliper to 1 from 2000.  Why was it 2000?
calada_scm <- lalonde_df_renamed %>%
  get_cal_matches(caliper = 0.5,
                  metric = METRIC,
                  cal_method = CAL_METHOD,
                  dist_scaling = DIST_SCALING,
                  est_method = "scm",
                  return = "sc_units",
                  #knn = 10,
                  num_bins = 5,
                  wider = FALSE)
calada_scm

get_att_point_est(calada_scm)

### Inference with the A-E method
calada_scm$hat_mu_0 <- 0

calada_scm_filtered <-
  calada_scm$result %>%
  filter(Z==F) %>%
  group_by(subclass) %>%
  filter(n() >= 2) %>%
  ungroup()

if ( FALSE ) {
  # Old bootstrap code--no longer in package
  sd_est = get_se_AE(calada_scm)
CI_lower[i] = mean_tilde_tau -1.96 *  sd_boot[i]
CI_upper[i] = mean_tilde_tau + 1.96 * sd_boot[i]
write_csv(boot_naive_res, file="sim_canonical_results/boot_naive_res.csv")
}

# Love plot vs. # co units added
set.seed(2)
love_plot(calada_scm,
          covs=paste0("X",5:8), B=NA) +
  scale_color_manual(
    values = wesanderson::wes_palette("Zissou1", 5)[c(5,3,2,1)]) +
  guides(color=F) +
  scale_x_continuous(breaks = c(174, 177, 180, 183))
ggsave("figures/lalonde_love.png", height=3, width=4)


if ( FALSE ) {
  # TODO: Add and export love_plot2 from package?  How is it different?

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
ggsave("figures/lalonde_love2.png", height=3, width=6)

}

# FSATT results -----------------------------------------------------------

feasible <- result_table( calada_scm, feasible_only = TRUE )
get_att_point_est(feasible)


# compare distances between units represented by average and sc weights
# scm_vs_avg_plot(feasible, DIST_SCALING, METRIC)
dist_density_plot(feasible, DIST_SCALING, METRIC) +
  theme(legend.position = "bottom")
ggsave("figures/lalonde_dist.png", height=3, width=4)

# compare ESS
ess_plot(feasible)
ggsave("figures/lalonde_ess.png", height=3, width=4)



# (optional) visualize density --------------------------------------------

require(latex2exp)
dist_matrix <- data.frame(t(as.matrix(attr(calada_scm, "dm_uncapped"))))
# rank each column
dm_col_sorted <- apply(dist_matrix, 2, sort)
dist_to_plot <- dm_col_sorted[1:3,] # choose top 3
tibble(d = as.numeric(as.matrix(dist_to_plot))) %>%
  filter(d < 1000) %>%
  ggplot(aes(d)) +
  geom_histogram(color="black", binwidth=0.3) +
  geom_vline(xintercept = 1, col="red") +
  theme_classic() +
  labs(y = "Count",
       x = TeX("$d(X_t, X_j)$"))
ggsave("figures/lalonde_calselect_top_3.png", height=2.5, width=5)






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
ggsave("figures/lalonde_att.png", height=3, width=6)
# [1] "FSATT: (339.098, 2851.353)"
# [1] "SATT: (76.115, 2612.28)"



