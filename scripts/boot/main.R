library(here)
load_all()
source(here("scripts/boot/boot_CSM_simulation_code.R"))

## Generate the bootstrap:
## Run the result
source(here("scripts/boot/boot_CSM_simulation_code.R"))
boot_otsu_wild <-
  boot_CSM(dgp_name="otsu",
           att0=T,
           I=100,
           B=100,
           mu_model="linear",
           boot_mtd="wild",
           n_split=2)

boot_otsu_A_E <-
  boot_CSM(dgp_name="otsu",
           att0=T,
           I=100,
           B=100,
           mu_model="linear",
           boot_mtd="A-E",
           n_split=2)
boot_otsu_naive_resid <-
  boot_CSM(dgp_name="otsu",
           att0=T,
           I=100,
           B=100,
           mu_model="linear",
           boot_mtd="naive-resid",
           n_split=2)
saveRDS(boot_otsu_wild, here("scripts/boot/output/boot_otsu_wild.rds"))
saveRDS(boot_otsu_A_E, here("scripts/boot/output/boot_otsu_A_E.rds"))
saveRDS(boot_otsu_naive_resid, here("scripts/boot/output/boot_otsu_naive_resid.rds"))


boot_otsu_naive_resid <- boot_otsu_naive_resid %>%
  rename(
    lower_wrong = lower,
    upper_wrong = upper,
    covered_wrong = covered
  )

# Recalculate lower and upper
boot_otsu_naive_resid <- boot_otsu_naive_resid %>%
  mutate(
    lower = lower_wrong / N_T,
    upper = upper_wrong / N_T
  )

# Recalculate covered
boot_otsu_naive_resid <- boot_otsu_naive_resid %>%
  mutate(
    covered = (lower < 0) & (0 < upper)
  )

# Load required libraries
library(here)
library(tidyverse)

# Read the saved results
boot_otsu_wild <- readRDS(here("scripts/boot/output/boot_otsu_wild.rds"))
boot_otsu_A_E <- readRDS(here("scripts/boot/output/boot_otsu_A_E.rds"))


create_bootstrap_comparison_plot(
  boot_otsu_wild = boot_otsu_wild,
  boot_otsu_A_E = boot_otsu_A_E,
  output_path =
    here("scripts/boot/figures/ci_comparison_plot.png"))



#### Experiments when there is overlap
## Run the result when overlap is low
source(here("scripts/boot/boot_CSM_simulation_code.R"))
boot_otsu_wild_low_overlap <-
  boot_CSM(dgp_name="otsu",
           att0=T,
           I=100,
           B=100,
           mu_model="linear",
           boot_mtd="wild",
           n_split=2,
           N1 = 20,
           N0 = 1000)

boot_otsu_A_E_low_overlap <-
  boot_CSM(dgp_name="otsu",
           att0=T,
           I=100,
           B=100,
           mu_model="linear",
           boot_mtd="A-E",
           n_split=2,
           N1 = 20,
           N0 = 1000)
saveRDS(boot_otsu_wild_low_overlap,
        here("scripts/boot/output/boot_otsu_wild_low_overlap.rds"))
saveRDS(boot_otsu_A_E_low_overlap,
        here("scripts/boot/output/boot_otsu_A_E_low_overlap.rds"))

boot_otsu_wild_low_overlap <- readRDS(here("scripts/boot/output/boot_otsu_wild_low_overlap.rds"))
boot_otsu_A_E_low_overlap <- readRDS(here("scripts/boot/output/boot_otsu_A_E_low_overlap.rds"))
source(here("scripts/boot/boot_CSM_simulation_code.R"))
create_bootstrap_comparison_plot(
    boot_otsu_wild = boot_otsu_wild_low_overlap,
    boot_otsu_A_E = boot_otsu_A_E_low_overlap,
    output_path =
      here("scripts/boot/figures/ci_comparison_plot_low_overlap.png"))



### Visulize realized DF ###
df_otsu <- gen_df_otsu(N = 100,K = 2)
ggplot(df_otsu, aes(x=X1,y=X2))+
  geom_point()
ggsave(here("scripts/boot/figures/otsu-rai-2d-x-plot.png"))

ggplot(df_otsu, aes(x=norm_X,y=m_X))+
  geom_point()
ggsave(here("scripts/boot/figures/otsu-rai-2d-y-plot.png"))


### For low overlap simulation ###
df_otsu <- gen_df_otsu(K = 2,N1 = 20,N0 = 1000)
ggplot(df_otsu, aes(x=norm_X,y=m_X))+
  geom_point()
