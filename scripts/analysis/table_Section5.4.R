

library(tidyr)
library(dplyr)
source( here::here( "scripts/analysis/MCSE_functions.R") )

# load data
R = 500
FNAME =
  here::here(
    paste0("data/outputs/A-E-overlap-by-prop-unif/",
           "A_E_toy_low_mid_high_R=",R,".csv")
  )
<<<<<<< HEAD
res <- read.csv(file = FNAME)
head( res )
=======
res <- read.csv(file =  FNAME) %>%
  rename(noise = "error") %>%
  mutate(error = noise + bias)
>>>>>>> 2d5c35c014983f579848ea864c575ec2742bbb0e

table_section_5_4 <- res %>%
  mutate(
    bias_and_var = (att_est - att_true),
    deg_overlap = factor(
      deg_overlap,
      # levels = c("low", "mid", "high")
      levels =  c("very_low", "low", "mid", "high", "very_high")
      )
  ) %>%
  group_by(deg_overlap) %>%
  summarise(
    SE_est = mean(se_AE),
    MCSE_SE_est = MCSE_bias(se_AE),
    # SE_True = mean(true_SE),
    # MCSE_SE_True = MCSE_bias(true_SE),
    SE_True = sd(noise),
    MCSE_SE_True = MCSE_SE(noise),
    # SE_samp_infl = sd(error + bias),
    # MCSE_SE_samp_infl = MCSE_SE(error + bias),
    SE_pop = sd(att_est),
    MCSE_SE_pop = MCSE_SE(att_est),
    N_C_tilde = N_C_tilde[1],
    mean_bias = mean(bias),
    MCSE_mean_bias = MCSE_bias(bias),
    coverage = mean(covered),
    coverage_with_true_SE = mean(covered_with_true_SE)
  ) %>%
  arrange(deg_overlap)

<<<<<<< HEAD
table_section_5_4
=======
write.csv(table_section_5_4,
          file = "tables/table_section_5_4.csv")
>>>>>>> 2d5c35c014983f579848ea864c575ec2742bbb0e
