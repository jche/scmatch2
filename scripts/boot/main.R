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
           boot_mtd="naive",
           n_split=2)
saveRDS(boot_otsu_wild, here("scripts/boot/output/boot_otsu_wild.rds"))
saveRDS(boot_otsu_A_E, here("scripts/boot/output/boot_otsu_A_E.rds"))
saveRDS(boot_otsu_naive_resid, here("scripts/boot/output/boot_otsu_naive_resid.rds"))


# Load required libraries
library(here)
library(tidyverse)
source(here("scripts/boot/boot_CSM_simulation_code.R"))
# Read the saved results
boot_otsu_wild <- readRDS(here("scripts/boot/output/boot_otsu_wild.rds"))
boot_otsu_A_E <- readRDS(here("scripts/boot/output/boot_otsu_A_E.rds"))

create_bootstrap_comparison_plot(
  boot_otsu_wild = boot_otsu_wild,
  boot_otsu_A_E = boot_otsu_A_E,
  output_path =
    here("scripts/boot/figures/ci_comparison_plot.png"))

##########
#### Experiments when debiasing is not used
source(here("scripts/boot/boot_CSM_simulation_code.R"))
boot_otsu_wild_undebiased <-
  boot_CSM(dgp_name="otsu",
           att0=T,
           I=100,
           B=100,
           boot_mtd="wild",
           debias = F)
source(here("scripts/boot/boot_CSM_simulation_code.R"))
boot_otsu_A_E_undebiased <-
  boot_CSM(dgp_name="otsu",
           att0=T,
           I=100,
           B=100,
           boot_mtd="A-E",
           debias = F)

saveRDS(boot_otsu_wild_undebiased, here("scripts/boot/output/boot_otsu_wild_undebiased.rds"))
saveRDS(boot_otsu_A_E_undebiased, here("scripts/boot/output/boot_otsu_A_E_undebiased.rds"))


########
#### Experiments when there is overlap -- not a focus now
## Run the result when overlap is low
source(here("scripts/boot/boot_CSM_simulation_code.R"))
boot_otsu_wild_low_overlap <-
  boot_CSM(dgp_name="otsu",
           att0=T,
           I=100,
           B=100,
           boot_mtd="wild",
           debias = F,
           N1 = 25,
           N0 = 1000)

boot_otsu_A_E_low_overlap <-
  boot_CSM(dgp_name="otsu",
           att0=T,
           I=100,
           B=100,
           boot_mtd="A-E",
           debias = F,
           N1 = 25,
           N0 = 1000)
saveRDS(boot_otsu_wild_low_overlap,
        here("scripts/boot/output/boot_otsu_wild_low_overlap.rds"))
saveRDS(boot_otsu_A_E_low_overlap,
        here("scripts/boot/output/boot_otsu_A_E_low_overlap.rds"))

source(here("scripts/boot/boot_CSM_simulation_code.R"))
boot_otsu_wild_low_overlap <- readRDS(here("scripts/boot/output/boot_otsu_wild_low_overlap.rds"))
boot_otsu_A_E_low_overlap <- readRDS(here("scripts/boot/output/boot_otsu_A_E_low_overlap.rds"))
create_bootstrap_comparison_plot(
    boot_otsu_wild = boot_otsu_wild_low_overlap,
    boot_otsu_A_E = boot_otsu_A_E_low_overlap,
    output_path =
      here("scripts/boot/figures/ci_comparison_plot_low_overlap.png"),
    show_paper_result = F)


########
## Run the bootstrap when block
########
source(here("scripts/boot/boot_CSM_simulation_code.R"))
boot_otsu_wild_moving_block <-
  boot_CSM(dgp_name="otsu",
           att0=T,
           I=100,
           B=100,
           boot_mtd="wild",
           debias = F,
           use_moving_block=T,
           block_size = 8)
mean(boot_otsu_wild_moving_block$covered)
saveRDS(boot_otsu_wild_moving_block,
        here("scripts/boot/output/boot_otsu_wild_moving_block.rds"))


### Visulize realized DF ###
df_otsu <- gen_df_otsu(N = 100,K = 2)
ggplot(df_otsu, aes(x=X1,y=X2))+
  geom_point(aes(colour = as.factor(Z) ))
ggsave(here("scripts/boot/figures/otsu-rai-2d-x-plot.png"))
## Next: visualize the overlapped KNNs

ggplot(df_otsu, aes(x=norm_X,y=m_X))+
  geom_point()
ggsave(here("scripts/boot/figures/otsu-rai-2d-y-plot.png"))


### For low overlap simulation ###
df_otsu <- gen_df_otsu(K = 2,N1 = 20,N0 = 1000)
ggplot(df_otsu, aes(x=norm_X,y=m_X))+
  geom_point()
