# scripts/toy_large_or_lm.R
#
# Generate a single large toy dataset (same DGP as sims-bias_mse/toy)
# and estimate the ATT using the full suite of methods.
#
# PURPOSE: Try to understand why the or_lm method had such good
# properties in the sim results sent to PA.
#
# CURRENT THEORY: It had the interaction of X1 and X2 in the outcome
# model.
#
# The DGP:
#   - X1, X2 ~ Uniform(0,1), controls clustered near center; treated near edges
#   - f0(X1,X2) = 20 * dmvnorm((X1,X2); mean=(0.5,0.5), Sigma=[[1,0.8],[0.8,1]])
#   - tau(X1,X2) = 3*X1 + 3*X2
#   - Y = f0(X1,X2) + Z*tau(X1,X2) + eps,  eps ~ N(0, f0_sd^2)

library(tidyverse)
library(mvtnorm)
library(optweight)
library(dbarts)
library(here)

devtools::load_all()
source(here::here("scripts/lib/wrappers.R"))
source(here::here("scripts/lib/sim_runner.R"))

set.seed(42)

# ------------------------------------------------------------------------------
# 1. Data generation
# ------------------------------------------------------------------------------

nc    <- 5000   # 10x the simulation default of 500
nt    <- 1000   # 10x the simulation default of 100
f0_sd <- 0.5
nbins <- 6      # same as simulation

df <- gen_df_adv(
  nc = nc,
  nt = nt,
  f0_sd = f0_sd,
  tx_effect_fun = function(X1, X2) { 3 * X1 + 3 * X2 },
  f0_fun = function(x, y) {
    matrix(c(x, y), ncol = 2) %>%
      dmvnorm(mean = c(0.5, 0.5),
              sigma = matrix(c(1, 0.8, 0.8, 1), 2)) * 20
  }
)

form <- as.formula("Z ~ X1 + X2")

cat(sprintf("Dataset: n=%d  (nt=%d, nc=%d)\n", nrow(df), sum(df$Z), sum(!df$Z)))

true_ATT <- df %>%
  filter(Z == 1) %>%
  summarize(att = mean(Y1 - Y0)) %>%
  pull(att)

cat(sprintf("True ATT: %.4f\n\n", true_ATT))

# ------------------------------------------------------------------------------
# 2. Scaling + feasibility filter
# ------------------------------------------------------------------------------

covs <- parse_form(form)$covs
dist_scaling <- df %>%
  summarize(across(all_of(covs), function(x) nbins / (max(x) - min(x))))

preds_csm <- get_cal_matches(data = df, treatment = "Z",
                             metric = "maximum", scaling = dist_scaling,
                             rad_method = "fixed", est_method = "average", k = 25)
preds_cem <- get_cem_matches(data = df, num_bins = nbins, est_method = "average",
                             return = "all")

ninf     <- length(attr(preds_csm, "unmatched_units"))
ninf_cem <- sum(df$Z) - length(attr(preds_cem, "feasible_units"))
df       <- df %>% filter(!id %in% attr(preds_csm, "unmatched_units"))

cat(sprintf("Infeasible units (CSM): %d,  (CEM): %d\n", ninf, ninf_cem))
cat(sprintf("Analysis n after filter: %d  (nt=%d)\n\n", nrow(df), sum(df$Z)))

# ------------------------------------------------------------------------------
# 3. Run all fast methods  (pass "all" for slow ones too)
# ------------------------------------------------------------------------------

results <- run_all_methods(
  df               = df,
  form             = form,
  dist_scaling     = dist_scaling,
  nbins            = nbins,
  selected_methods = "fast"
) %>%
  mutate(
    true_ATT = true_ATT,
    bias     = ATT_est - true_ATT
  )

cat("\n--- Results ---\n")
cat(sprintf("True ATT: %.4f\n\n", true_ATT))
print( results %>%
        dplyr::select(method, ATT_est, bias, secs), n = Inf)



# Look at linear regression on this dataset ----

df
#install.packages("hexbin")
library( hexbin )
ggplot( df, aes( X1, X2 ) ) +
  geom_hex( bins = 30 ) +
  coord_fixed()

df0 = filter( df, Z == 0 )
df0$X2ct = cut( df0$X2, breaks = 10 )
ggplot( df0, aes( X1, Y ) ) +
  facet_wrap( ~ X2ct ) +
  geom_point()
ggplot( df0, aes( X2, Y ) ) +
  geom_point()


M = lm( Y ~ (X1 + X2)^2, data = filter( df, Z == 0 ) )
library( broom )
tidy( M )
df1 = dplyr::filter( df, Z == 1 )
df1$pdY0 = predict( M, newdata=df1 )
mean( df1$Y - df1$pdY0 )

get_att_or_lm( df, form=Y ~ X1+X2 )
get_att_or_lm( df, form=Y ~ (X1+X2)^2 )



