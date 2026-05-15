
# Code to explore single datasets and see how the functions all work


suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(here)
})
source(here::here("scripts/sims-variance-multi/0_utils.R"))


# Make data ----
if ( FALSE ) {
  debugonce( make_df_multi )
}

df <- make_df_multi( error_type = "het",
                     sigma1_extra = 0,
                     prop_nc_unif = 0.05 )

head( df )
stats <- df %>%
  mutate( err0 = Y0 - Y0_denoised,
          err1 = Y1 - Y1_denoised ) %>%
  summarize( sdY0 = sd( Y0 ),
             sd0 = sd( err0 ),
             sd1 = sd( err1 ),
             mn_sigma0 = sqrt( mean( sigma0^2 ) ) )
stats




# Match data ---
mtch <-  get_cal_matches(
  data       = df,
  Z ~ X1 + X2,
  rad_method = "adaptive",
  caliper = 0.1,
  scaling    = 1,
  k          = 1,
  warn       = FALSE,
  est_method = "scm"
)
mtch

summary( mtch )

true_satt <- df %>%
  filter(Z == 1) %>%
  summarise(att = mean(Y1 - Y0)) %>%
  pull(att)
true_satt

vrs <- get_finite_variance( mtch,
                     use_common_variance = TRUE,
                     homoskedastic = FALSE )
vrs


# Compare estimated noise to the true noise
bind_cols( vrs, stats ) %>%
  mutate( S1 = sqrt(S1_sq),
          S0 = sqrt(S0_sq) ) %>%
  dplyr::select( -SE, -V_E, -cov_w_s, -sigma_hat, -S1_sq, -S0_sq )



