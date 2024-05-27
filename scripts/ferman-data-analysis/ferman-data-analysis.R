library(haven)

data_ferman <- read_dta(
  file = "./data/inputs/Final.dta"
)
colnames(data_ferman)
hist(data_ferman$y2008 - data_ferman$OBJETIVA2008)
ferman_for_analysis <-
  data_ferman %>%
  select(starts_with("y20"),
         Control,
         UF) %>%
  filter(if_all(starts_with("y20"), ~ !is.na(.x))) %>%
  filter(!is.na(Control)) %>%
  mutate(Y = y2010,
         Z = Control)

test_UF <- ferman_for_analysis %>%
  group_by(Control, UF) %>%
  summarise(n())

hist(ferman_for_analysis$y2007)
# First, plot a histogram of 1nn  to select scaling and caliper
CALIPER <- 0.02
METRIC <- "maximum"   # "maximum", "euclidean", "manhattan"
CAL_METHOD <- "ada"   # "ada", "cem"
DIST_SCALING <- 0.3
ferman_scm <- ferman_for_analysis %>%
  get_cal_matches(
    covs = c("y2007", "y2008", "y2009"),
    treatment = "Z",
    caliper = CALIPER,
                metric = METRIC,
                cal_method = CAL_METHOD,
                dist_scaling = DIST_SCALING,
                est_method = "scm",
                return = "sc_units",
                knn = 10,         # for ada
                num_bins = 5,     # for cem
                wider = F)
get_att_ests(ferman_scm)

dist_matrix <- data.frame(t(as.matrix(attr(ferman_scm, "dm_uncapped"))))
# rank each column
dm_col_sorted <- apply(dist_matrix, 2, sort)
dist_to_plot <- dm_col_sorted[1:3,] # choose top 3
tibble(d = as.numeric(as.matrix(dist_to_plot))) %>%
  filter(d < 1000) %>%
  ggplot(aes(d)) +
  geom_histogram(color="black") +
  geom_vline(xintercept = 0.02, col="red") +
  theme_classic() +
  labs(y = "Count",
       x = TeX("$d(X_t, X_j)$"))


ferman_dm <-
  gen_dm(df = ferman_for_analysis,
       covs = c("y2007", "y2008", "y2009"),
       treatment = "Control",
       scaling = 0.3,
       metric = "maximum"
)

ferman_dm[1:5, 1:5]

res <- fit_CSM(df = ferman_for_analysis,
               covs = c("y2007", "y2008", "y2009"),
               treatment = "Control",
               metric = "maximum",
               caliper = 1,
               rad_method = "adaptive",
               est_method = "scm",
               return = "agg_co_units",
               dist_scaling = 1)

