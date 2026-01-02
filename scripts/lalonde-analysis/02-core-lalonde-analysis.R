# scripts/lalonde-analysis/02-core-lalonde-analysis.R
# Conduct the core analysis for prepping summary results

library(CSM)
library(tidyverse)
library(here)

lalonde_df <- readRDS(here::here("scripts/lalonde-analysis/data/lalonde_for_analysis.rds"))
lalonde_df$HS = ifelse( lalonde_df$X6 >= 12, 1, 0 )
lalonde_df$COLL = ifelse( lalonde_df$X6 >= 16, 1, 0 )

# Set X1 to X4 to numeric
lalonde_df <- lalonde_df %>%
  mutate(
    X1 = as.numeric(X1),
    X2 = as.numeric(X2),
    X3 = as.numeric(X3),
    X4 = as.numeric(X4)
  )

dim( lalonde_df )


skimr::skim( lalonde_df )

# EXPLORING WHICH VARIABLES ARE IMPORTANT ----

library(broom)
library(vip)
library(twang)

set.seed(232342)

ps_fit <- glm(
  Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8,
  data = lalonde_df,
  family = binomial()
)

# propensity scores
lalonde_df <- lalonde_df %>%
  mutate(pscore_glm = predict(ps_fit, type = "response"))

lalonde_df %>%
  group_by( Z ) %>%
  summarise( q5 = quantile( pscore_glm, 0.05 ),
             q25 = quantile( pscore_glm, 0.25 ),
             q50 = quantile( pscore_glm, 0.5 ),
             q75 = quantile( pscore_glm, 0.75 ),
             q95 = quantile( pscore_glm, 0.95 ),
             max = max( pscore_glm ),
  )
mean( lalonde_df$pscore_glm[ lalonde_df$Z == 0 ] <= 0.0124 )


# importance by absolute standardized coefficient (|beta| on z-scored predictors)
X_mat <- model.matrix(ps_fit)[, -1, drop = FALSE]
beta  <- coef(ps_fit)[-1]
imp   <- tibble(term = colnames(X_mat), beta = beta) %>%
  mutate(
    x_sd = apply(X_mat, 2, sd),
    beta_std = beta * x_sd,
    importance = abs(beta_std)
  ) %>%
  arrange(desc(importance))
imp

# quick plot
vip::vip(ps_fit)

round( summary( lalonde_df$pscore_glm ), digits = 4 )

lalonde_df$logit_p = predict( ps_fit, type="link" )
0.2 * sd( lalonde_df$logit_p[ lalonde_df$Z==1] )

ggplot( lalonde_df, aes( x=logit_p, fill=as.factor(Z) ) ) +
  geom_boxplot( alpha=0.5 ) +
  theme_minimal() +
  labs( fill = "Treatment" ) +
  ggtitle( "Propensity Score Distribution by Treatment Group" )


# Trying same thing with twang ----

if ( FALSE ) {

  # twang GBM propensity score model
  ps_twang <- ps(
    formula = Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8,
    data = as.data.frame( lalonde_df ),
    n.trees = 5000,
    interaction.depth = 3,
    shrinkage = 0.01,
    perm.test.iters = 0,
    stop.method = c("es.mean", "ks.max"),
    estimand = "ATE",
    verbose = FALSE
  )

  # pick stopping rule (use es.mean by default)
  stop_rule <- "es.mean"

  # get propensity scores
  lalonde_df$pscore <- twang::get.weights(ps_twang, stop.method = stop_rule)

  # covariate importance (relative influence from GBM)
  var_imp <- summary(ps_twang, stop.method = stop_rule)$rel.inf %>%
    as_tibble() %>%
    rename(var = var, rel_inf = rel.inf) %>%
    arrange(desc(rel_inf))

  var_imp

  # optional: plot balance + influence
  plot(ps_twang, plots = 1, stop.method = stop_rule)  # balance
  plot(ps_twang, plots = 3, stop.method = stop_rule)  # relative influence

}



# MATCHING SETTINGS

METRIC <- "maximum" # "euclidean"
CAL_METHOD <- "adaptive" # The paper uses adaptive caliper: c_t = max{1, d_t}


# Scaling setup matches the paper's V matrix:
# V = diag {K, K, K, K, 1/3, 1, 1/5000, 1/5000}, where K=1000.
DIST_SCALING <- tibble(
  X1 = 1000, # % Black (Exact Match)
  X2 = 1000, # % Hispanic (Exact Match)
  X3 = 1000, # % Married (Exact Match)
  X4 = 1000, # % No degree (Exact Match)
  X5 = 1/3,  # Average age (Caliper of 3 years)
  X6 = 1/2,    # Average years of education (Caliper of 2 years)
  X7 = 1/5000, # Average earnings, 1974 (Caliper of $5,000)
  X8 = 1/1400,  # Average earnings, 1975 (Caliper of $5,000)
  COLL = 5/7,
  logit_p = 1 / 0.35
)


# Run matching
lalonde_scm <- lalonde_df %>%
  get_cal_matches(
    form = Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + COLL + logit_p,
    caliper = 1,
    metric = METRIC,
    rad_method = CAL_METHOD,
    scaling = DIST_SCALING,
    est_method = "scm"
  )

lalonde_scm

cdp <- caliper_distance_plot( lalonde_scm, target_percentile = 2/3 )
cdp



# based on plot
CALIPER = 0.552

lalonde_scm <- update_matches( lalonde_scm, lalonde_df, caliper = CALIPER )
lalonde_scm

