# scripts/ferman-analysis/ferman-data-analysis.R
library(haven)
library(CSM)
library(latex2exp)
library(tidyverse)

options(list(dplyr.summarise.inform = FALSE))
theme_set( theme_classic() )

# Load the data
ferman_for_analysis <-
  readRDS( here::here( "scripts/ferman-analysis/data/ferman_for_analysis.rds" ) )
head( ferman_for_analysis )

# Drop unneeded variables
ferman_for_analysis$UF = NULL
ferman_for_analysis$y2011 = NULL
ferman_for_analysis$y2012 = NULL
ferman_for_analysis$Control = NULL
ferman_for_analysis$Y = NULL

names( ferman_for_analysis )

# The main matching analysis ----
c <- 0.35
covariate_caliper <- c(rep(0.2, 3), 1/1000)
ferman_scm <- get_cal_matches( ferman_for_analysis,
    covs = c("y2007", "y2008", "y2009", "is_sao_paolo"),
    treatment = "Z",
    caliper = c,
    metric = "maximum",   # "maximum", "euclidean", "manhattan"
    rad_method = "adaptive",
    scaling = 1/covariate_caliper,
    est_method = "scm" )
ferman_scm

result_table(ferman_scm, nonzero_weight_only = TRUE ) %>%
  rename( isp = is_sao_paolo )

estimate_ATT( ferman_scm, "Z", "y2010" )


## Number of used controls

summary( ferman_scm )



#####
# Get the balance table ----
#####
## TODO: fix this piece
# d <- feasible
# metric <- "maximum"
# sc_dists <-
#   gen_dm(df = d %>%
#            CSM:::agg_sc_units(),
#          covs = covs,
#          treatment = "Z",
#          scaling = scaling,
#          metric = metric)%>%
#   diag()
#
#
#
# avg_dists <-
#   gen_dm(df = d %>%
#            CSM:::agg_avg_units(),
#          covs = covs,
#          treatment = "Z",
#          scaling = scaling,
#          metric = metric) %>%
#   diag()
#
# nn_dists <-
#   gen_dm(df = d %>%
#            group_by(Z,subclass) %>%
#            filter(dist == min(dist)) %>%
#            slice(1) %>%
#            mutate(weights = 1) %>%
#            ungroup() %>%
#            CSM:::agg_avg_units(),
#          covs = covs,
#          treatment = "Z",
#          scaling = scaling,
#          metric = metric) %>%
#   diag()
# mean(sc_dists)
# mean(avg_dists)
# mean(nn_dists)

# res_list <-
#   list(sc = sc_dists,
#        avg = avg_dists,
#        nn = nn_dists)
# print(tibble(
#   method = c("SCM", "Average", "1-NN"),
#   mean = map_dbl(res_list, ~round(mean(.), 3)),
#   median = map_dbl(res_list, ~round(median(.), 3))
# ))


######
# Histogram ----
######
feasible_subclasses <- feasible_unit_subclass(ferman_scm)
(n_feasible <- length(feasible_subclasses) )

matched_controls <- result_table( ferman_scm, "all" ) %>%
  group_by(subclass) %>%
  summarise(n_controls = n()-1) %>%
  filter(subclass %in% feasible_subclasses)

sum(matched_controls$n_controls)
hist_matched_controls <-
  ggplot(matched_controls, aes(x = n_controls)) +
  geom_histogram(binwidth = 3, fill = "skyblue", color = "white") +
  labs(
    x = "Number of Matched Controls",
    y = "Frequency"
  ) +
  theme_minimal()
hist_matched_controls

matched_controls_nonzero_weights <-
  ferman_scm$matches %>%
  bind_rows() %>%
  filter(weights > 0) %>%
  group_by(subclass) %>%
  summarise(n_controls = n()-1) %>%
  filter(subclass %in% feasible_subclasses)
sum(matched_controls_nonzero_weights$n_controls)

hist_matched_controls_nonzero_weights<-
  ggplot(matched_controls_nonzero_weights, aes(x = n_controls)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "white") +
  labs(
    x = "Number of Matched Controls",
    y = "Frequency"
  ) +
  theme_minimal()
p <- gridExtra::grid.arrange(
  hist_matched_controls,
  hist_matched_controls_nonzero_weights,
  ncol=2
)
ggsave(
  here::here( "figures/hist-n-co.png" ),
  plot=p,
  width=8.1,
  height=5.3
)


#####
# SATT plot ----
######

ggd_att <-
  result_table( ferman_scm ) %>%
  left_join( caliper_table( ferman_scm ),
             by = c( "id", "subclass" ) ) %>%
  group_by(subclass) %>%
  summarize(adacal = last(na.omit(adacal)),
            tx = y2010[2] - y2010[1]) %>%
  arrange(adacal) %>%
  mutate(order = 1:n(),
         cum_avg = cumsum(tx) / order ) %>%
  slice((n_feasible):n())
ggd_att

plot_max_caliper_size <-
  ggd_att %>%
  ggplot(aes(x=order, y=adacal)) +
  geom_line(alpha=0.5) +
  geom_point(size=3) +
  theme_classic() +
  labs(y = "Maximum caliper size used",
       x = "Total number of treated units used") +
  expand_limits(color=1)

feasible_w_adacal <-
  result_table(ferman_scm) %>%
  left_join( caliper_table( ferman_scm ),
             by = c( "id", "subclass" ) ) %>%
  group_by(subclass)

# feasible_w_adacal <-
#   attr(ferman_scm, "scweights") %>%
#   bind_rows() %>%
#   left_join(attr(ferman_scm, "adacalipers"),
#             by="id")

# make all subclass in tmp to have the same adacal
# currently some subclass has NA in adacal
feasible_w_adacal <-
  feasible_w_adacal %>%
  group_by(subclass) %>%
  mutate(adacal =
           ifelse(is.na(adacal),
                  first(na.omit(adacal)),
                  adacal)) %>%
  ungroup() %>%
  arrange(desc(adacal))


feasible_w_adacal$hat_mu_0 <- 0
# Number of unique subclass levels
n_unique_subclass <-
  length(unique(feasible_w_adacal$subclass))

# Pre-allocate a numeric vector to store the results
se_AEs <- numeric(n_unique_subclass)

# Iterate over the range from the number of unique subclass levels down to a certain number (e.g., 45)
for (i in n_unique_subclass:n_feasible) {
  # Filter the dataframe to only include rows where the subclass value is within the top 'i' lowest adacal values
  top_subclasses <- feasible_w_adacal %>%
    group_by(subclass) %>%
    summarise(min_adacal = min(adacal)) %>%
    arrange(min_adacal) %>%
    slice(1:i) %>%
    pull(subclass)

  df_curr <- feasible_w_adacal %>%
    filter(subclass %in% top_subclasses)

  # Get SE
  source(here::here("R/estimate.R"))

  preds_csm_filtered <- df_curr %>%
    filter(Z == 0) %>%
    group_by(subclass) %>%
    filter(n() >= 2) %>%
    ungroup()

  # 3. For each filtered subclass,
  #   calculate the cluster residual s_j^2 = se(debiased_units)
  weighted_var <- function(x, wt) {
    n <- length(x)
    wm <- weighted.mean(x, wt)
    sum(wt * (x - wm)^2) * n / (n-1)
  }

  weighted_se <- function(x, wt) {
    sqrt(weighted_var(x, wt) / length(x))
  }

  cluster_var_df <-
    preds_csm_filtered %>%
    group_by(subclass) %>%
    summarise(nj = n(),
              w_nj = ess(weights),
              var_cluster = var(y2010), .groups="drop")

  weighted_var_df <- cluster_var_df %>%
    summarise(weighted_var = weighted.mean(var_cluster, w = w_nj), .groups="drop")
  sigma_hat <- sqrt(weighted_var_df$weighted_var)

  # 5. calculate N_T and N_C
  Ns <- calc_N_T_N_C(df_curr, treatment="Z")

  # 6. Calculate the variance of the estimator
  sd_curr <- sqrt((1/Ns$N_T + 1/Ns$N_C_tilde)) * sigma_hat

  se_AEs[i] <- sd_curr
}

# Output the results
foo <-
  ggd_att %>%
  ggplot(aes(x=order, y=cum_avg)) +
  geom_line(alpha=0.5) +
  geom_point(size=3) +
  theme_classic() +
  labs(y = "Cumulative ATT Estimate",
       x = "Total number of treated units used",
       color = "Maximum caliper size used") +
  expand_limits(color=1)

p <- foo +
  geom_hline(yintercept=0, lty="dotted") +
  theme(legend.direction="horizontal",
        legend.position = c(0.5, 0.85),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  labs(y = "Cumulative ATT Estimate",
       x = "Total number of treated units used",
       color = "Maximum \ncaliper \nsize used    ")+
  ylim(c(-0.2,0.5))
foo$data$cum_avg
se_AEs
plot_SATT <- p+
  geom_errorbar(
    aes(ymin=cum_avg-1.96*se_AEs[n_feasible:54],
        ymax=cum_avg+1.96*se_AEs[n_feasible:54]),
    width = 0.5,
    linewidth=1
  ) +
  ylim(c(-0.1,0.2))


plot_max_caliper_size
plot_SATT



# Make a combined plot for the paper
plot_all <- gridExtra::grid.arrange(
  plot_max_caliper_size,
  plot_SATT,
  ncol=2
)

ggsave(
  filename="figures/fsatt-ferman.png",
  plot=plot_all,
  width=9.7,
  height=4.5
)




#########
# Estimate using average and 1-nn with MatchIt ----
#########

library(MatchIt)
covs <- c("y2007", "y2008", "y2009", "is_sao_paolo")
covs_backtick <-
  sapply(covs, function(x) paste0("`", x, "`"))
formula <-
  as.formula(
    paste("Z ~",
          paste(covs_backtick, collapse = " + ")))
ferman_1nn_match <- matchit(
  formula,
  data = ferman_for_analysis,
  method = "nearest", # used to be 'average'?
  distance = "euclidean")

summary(ferman_1nn_match)
m.data <- match.data(ferman_1nn_match)

#Linear model with covariates
formula_outcome <-
  as.formula(
    paste("y2010 ~ Z*(",
          paste(covs_backtick, collapse = " + "),
          ")"
        )
    )
fit1 <- lm(formula_outcome,
           data = m.data,
           weights = weights)
marginaleffects::avg_comparisons(
  fit1, variables = "Z",
        vcov = ~subclass,
        newdata = subset(m.data, Z == 1),
        wts = "weights")


## Average
ferman_scm_avg <- ferman_for_analysis %>%
  get_cal_matches(
    covs = covs,
    treatment = "Z",
    caliper = c,
    metric = "maximum",   # "maximum", "euclidean", "manhattan"
    rad_method = "adaptive",
    scaling = 1/covariate_caliper,
    est_method = "average" )

feasible_w_adacal <-
  result_table(ferman_scm_avg) %>%
  left_join( caliper_table( ferman_scm_avg ),
             by = c( "id", "subclass" ) ) %>%
  group_by(subclass)


feasible_w_adacal <-
  feasible_w_adacal %>%
  group_by(subclass) %>%
  mutate(adacal =
           ifelse(is.na(adacal),
                  first(na.omit(adacal)),
                  adacal)) %>%
  ungroup() %>%
  arrange(desc(adacal))


feasible_w_adacal$hat_mu_0 <- 0
# Number of unique subclass levels
n_unique_subclass <-
  length(unique(feasible_w_adacal$subclass))

# Pre-allocate a numeric vector to store the results
se_AEs <- numeric(n_unique_subclass)


  df_curr <- feasible_w_adacal

  # Get SE
  source(here::here("R/estimate.R"))

  preds_csm_filtered <- df_curr %>%
    filter(Z == 0) %>%
    group_by(subclass) %>%
    filter(n() >= 2) %>%
    ungroup()

  # 3. For each filtered subclass,
  #   calculate the cluster residual s_j^2 = se(debiased_units)
  weighted_var <- function(x, wt) {
    n <- length(x)
    wm <- weighted.mean(x, wt)
    sum(wt * (x - wm)^2) * n / (n-1)
  }

  weighted_se <- function(x, wt) {
    sqrt(weighted_var(x, wt) / length(x))
  }

  cluster_var_df <-
    preds_csm_filtered %>%
    group_by(subclass) %>%
    summarise(nj = n(),
              w_nj = ess(weights),
              var_cluster = var(y2010), .groups="drop")

  weighted_var_df <- cluster_var_df %>%
    summarise(weighted_var = weighted.mean(var_cluster, w = w_nj), .groups="drop")
  sigma_hat <- sqrt(weighted_var_df$weighted_var)

  # 5. calculate N_T and N_C
  Ns <- calc_N_T_N_C(df_curr, treatment="Z")

  # 6. Calculate the variance of the estimator
  sd_curr <- sqrt((1/Ns$N_T + 1/Ns$N_C_tilde)) * sigma_hat

  (se_AE_avg <- sd_curr)
  get_att_point_est(df_curr, outcome = "y2010")

