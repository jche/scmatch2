library(here)
load_all()
source(here("scripts/boot/boot_CSM_simulation_code.R"))

### Visulize realized DF
df_otsu <- gen_df_otsu(N = 100,K = 2)
ggplot(df_otsu, aes(x=X1,y=X2))+
  geom_point()
ggsave(here("scripts/boot/figures/otsu-rai-2d-x-plot.png"))

ggplot(df_otsu, aes(x=norm_X,y=m_X))+
  geom_point()
ggsave(here("scripts/boot/figures/otsu-rai-2d-y-plot.png"))

## Generate the bootstrap:
# The main procedure:
#   Repeat I times of DGP --> match --> getting bootstrapped intervals

# Find the code of the bootstrap procedure that calls generate_one_otsu
# Result: get_df_scaling_from_dgp_name in here("scripts/boot/boot_CSM_simulation_code.R")

boot_otsu_bayesian <-
  boot_CSM(dgp_name="otsu",
           att0=T,
           I=100,
           B=100,
           mu_model="linear",
           boot_mtd="Bayesian",
           n_split=2)
mean(boot_otsu_bayesian$covered)


boot_otsu_A_E <-
  boot_CSM(dgp_name="otsu",
           att0=T,
           I=100,
           B=100,
           mu_model="linear",
           boot_mtd="A-E",
           n_split=2)
mean(boot_otsu_A_E$covered)


test_boot_by_resids()
test_boot_bayesian()

