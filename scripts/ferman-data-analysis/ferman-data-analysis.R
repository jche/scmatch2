library(haven)
library(CSM)
library(latex2exp)

data_ferman <- read_dta( file = here::here( "data/inputs/Final.dta" ) )

colnames(data_ferman)
hist(data_ferman$y2008 -
       data_ferman$OBJETIVA2008)
# The above shows that y2008 and OBJETIVA2008 are the same variable
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
ferman_scm <- ferman_for_analysis %>%
  get_cal_matches(
    covs = c("y2007", "y2008", "y2009"),
    treatment = "Z",
    caliper = 0.04,
    metric = "maximum",   # "maximum", "euclidean", "manhattan"
    rad_method = "adaptive",
    dist_scaling = 1,
    est_method = "scm",
    return = "sc_units",
    knn = 10,         # for ada
    num_bins = 5,     # for cem
    wider = F)
dist_matrix <-
  data.frame(
    t(as.matrix(attr(ferman_scm, "dm_uncapped"))))
# rank each column
dm_col_sorted <- apply(dist_matrix, 2, sort)
dist_to_plot <- dm_col_sorted[1:3,] # choose top 3
tibble(d = as.numeric(as.matrix(dist_to_plot))) %>%
  filter(d < 1000) %>%
  ggplot(aes(d)) +
  geom_histogram(color="black", binwidth=0.01) +
  geom_vline(xintercept = 0.04, col="red") +
  theme_classic() +
  labs(y = "Count",
       x = TeX("$d(X_t, X_j)$"))


####
### Love plot
scmatches <- ferman_for_analysis %>%
  CSM::gen_matches(
    covs = c("y2007", "y2008", "y2009"),
    treatment = "Z",
    scaling = 1,
    metric = "maximum",
    caliper = 0.05,
    rad_method = "adaptive")
matched_co_unit1 <-  scmatches$matches[[1]]

scweights <-
  est_weights(
    df = ferman_for_analysis,
    covs=c("y2007", "y2008", "y2009"),
    matched_gps = scmatches$matches,
    dist_scaling = 1,
    est_method = "scm",
    metric = "maximum")
matched_co_unit1_with_weights <-
  scweights[[1]]
ferman_scm$weights


ferman_scm <- ferman_for_analysis %>%
  get_cal_matches(
    covs = c("y2007", "y2008", "y2009"),
    treatment = "Z",
    caliper = 0.04,
    metric = "maximum",   # "maximum", "euclidean", "manhattan"
    rad_method = "adaptive",
    dist_scaling = 1,
    est_method = "scm",
    return = "sc_units",
    knn = 10,         # for ada
    num_bins = 5,     # for cem
    wider = F)
# source("R/diagnostic_plots.R")
love_plot(res = ferman_scm,
           covs = c("y2007", "y2008", "y2009"),
           B=NA)

## ESS plots
scweights_df <-
  attr(ferman_scm, "scweights") %>%
  bind_rows()
subclass_feasible <-
  scweights_df %>%
  filter(id %in%
           attr(ferman_scm, "feasible_units") )%>%
  pull(subclass)

feasible <-
  scweights_df %>%
  filter(subclass %in%
           subclass_feasible)
ess_plot(feasible)

#####
## Get the balance table

# compare distances between units represented by average and sc weights
d <- feasible
metric <- "maximum"
dist_scaling <- 1
source("R/distance.R")
source("R/estimate.R")

sc_dists <-
  gen_dm(df = d %>%
           agg_sc_units(),
         covs = c("y2007", "y2008", "y2009"),
         treatment = "Z",
         scaling = dist_scaling,
         metric = metric)%>%
  diag()

avg_dists <-
  gen_dm(df = d %>%
           agg_avg_units(),
         covs = c("y2007", "y2008", "y2009"),
         treatment = "Z",
         scaling = dist_scaling,
         metric = metric) %>%
  diag()

nn_dists <-
  gen_dm(df = d %>%
           group_by(Z,subclass) %>%
           filter(dist == min(dist)) %>%
           slice(1) %>%
           mutate(weights = 1) %>%
           ungroup() %>%
           agg_avg_units(),
         covs = c("y2007", "y2008", "y2009"),
         treatment = "Z",
         scaling = dist_scaling,
         metric = metric) %>%
  diag()


mean(sc_dists)
mean(avg_dists)
mean(nn_dists)


res_list <-
  list(sc = sc_dists,
       avg = avg_dists,
       nn = nn_dists)
print(tibble(
  method = c("SCM", "Average", "1-NN"),
  mean = map_dbl(res_list, ~round(mean(.), 3)),
  median = map_dbl(res_list, ~round(median(.), 3))
))


#####
## SATT plot
######
ferman_scm$id <-
  gsub("*_syn","",ferman_scm$id)
ferman_scm$id <- as.integer(ferman_scm$id)
foo <- satt_plot4(ferman_scm, B=NA)

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
source("./R/bootstrap.R")

feasible_subclasses <- attr(ferman_scm, "feasible_subclasses")
n_feasible <- length(feasible_subclasses)
feasible_w_adacal <- attr(ferman_scm, "scweights") %>%
  bind_rows() %>%
  left_join(attr(ferman_scm, "adacalipers"), by="id")
# make all subclass in tmp to have the same adacal
# currently some subclass has NA in adacal
feasible_w_adacal <- feasible_w_adacal %>%
  group_by(subclass) %>%
  mutate(adacal = ifelse(is.na(adacal), first(na.omit(adacal)), adacal)) %>%
  ungroup() %>%
  arrange(desc(adacal))


feasible_w_adacal$hat_mu_0 <- 0
# Number of unique subclass levels
n_unique_subclass <- length(unique(feasible_w_adacal$subclass))

# Pre-allocate a numeric vector to store the results
se_AEs <- numeric(n_unique_subclass)

# Iterate over the range from the number of unique subclass levels down to a certain number (e.g., 45)
for (i in n_unique_subclass:45) {
  # Filter the dataframe to only include rows where the subclass value is within the top 'i' lowest adacal values
  top_subclasses <- feasible_w_adacal %>%
    group_by(subclass) %>%
    summarise(min_adacal = min(adacal)) %>%
    arrange(min_adacal) %>%
    slice(1:i) %>%
    pull(subclass)

  df_curr <- feasible_w_adacal %>%
    filter(subclass %in% top_subclasses)

  # Calculate the statistic and save it to the vector
  se_AEs[i] <- get_se_AE(df_curr)
}

# Output the results
se_AEs
p+
  geom_errorbar(aes(ymin=cum_avg-1.96*se_AEs[45:54],
                    ymax=cum_avg+1.96*se_AEs[45:54]),
                width = 0.5,
                linewidth=1) +
  geom_point(aes(color=adacal),
             size=3)

p+
  geom_errorbar(aes(ymin=cum_avg-1.96*se_AEs[54],
                    ymax=cum_avg+1.96*se_AEs[54]),
                width = 0.5,
                linewidth=1) +
  geom_point(aes(color=adacal),
             size=3)

#####
# DRAFT
####
diff_scm_co_and_tx<-
  get_diff_scm_co_and_tx(
    res = ferman_scm,
    covs = c("y2007", "y2008", "y2009"))
love_plot_df <-
  create_love_plot_df(res = ferman_scm,
                      covs = c("y2007", "y2008", "y2009"))


dist_matrix <- data.frame(t(as.matrix(attr(ferman_scm, "dm_uncapped"))))

get_att_ests(ferman_scm)
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

