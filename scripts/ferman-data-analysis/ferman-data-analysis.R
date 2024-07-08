library(haven)
library(CSM)
library(latex2exp)

data_ferman <- read_dta(
  file = "./data/inputs/Final.dta"
)
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
         Z = Control,
         is_sao_paolo = UF == 35)%>%
  mutate(`X1 (2007 score)` = y2007,
         `X2 (2008 score)` = y2008,
         `X3 (2009 score)` = y2009,
         `X4 (pct Sao Paolo)` = is_sao_paolo)
# UF == 35 means Sao Paolo
# UF == 33 means Rio


proportion_by_UF <- ferman_for_analysis %>%
  group_by(Control, is_sao_paolo) %>%
  summarise(n())

# First, plot a histogram of 1nn  to select scaling and caliper
c <- 0.35
covariate_caliper <- c(rep(0.2, 3), 1/1000)
ferman_scm <- ferman_for_analysis %>%
  get_cal_matches(
    covs = c("y2007", "y2008", "y2009", "is_sao_paolo"),
    treatment = "Z",
    caliper = c,
    metric = "maximum",   # "maximum", "euclidean", "manhattan"
    rad_method = "adaptive",
    dist_scaling = 1/covariate_caliper,
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
plot_dm <- function(dist_to_plot){
  tibble(d = as.numeric(as.matrix(dist_to_plot))) %>%
    ggplot(aes(d)) +
    geom_histogram(color="black", binwidth=0.03) +
    geom_vline(xintercept = c, col="red") +
    theme_classic() +
    labs(y = "Count",
         x = TeX("$d(X_t, X_j)$")) +
    xlim(c(0,1))
}
dist_to_plot <- dm_col_sorted[1:3,] # choose top 3
plot_dm(dist_to_plot)
# Plot top 1, 2, 3 respectively
library(gridExtra)
dm_top_1 <- dm_col_sorted[1,]
dm_top_2 <- dm_col_sorted[2,]
dm_top_3 <- dm_col_sorted[3,]
plot_dm_all <- gridExtra::grid.arrange(
  plot_dm(dm_top_1),
  plot_dm(dm_top_2),
  plot_dm(dm_top_3)
)

ggsave(
  filename="writeup/figures/hist-top-k-distances.png",
  plot=plot_dm_all,
  width=5.3,
  height=7.7
)

####
### Love plot
scaling <- 1/covariate_caliper
covs <-c("X1 (2007 score)", "X2 (2008 score)",
"X3 (2009 score)", "X4 (pct Sao Paolo)")
# covs <- starts_with("X")
# covs <- c("X1", "X2", "X3", "X4")
# scmatches <- ferman_for_analysis %>%
#   CSM::gen_matches(
#     covs = covs,
#     treatment = "Z",
#     scaling = scaling,
#     metric = "maximum",
#     caliper = 0.015,
#     rad_method = "adaptive")
# matched_co_unit1 <-  scmatches$matches[[1]]
#
# scweights <-
#   est_weights(
#     df = ferman_for_analysis,
#     covs= covs,
#     matched_gps = scmatches$matches,
#     dist_scaling = scaling,
#     est_method = "scm",
#     metric = "maximum")
# matched_co_unit1_with_weights <-
#   scweights[[1]]
# ferman_scm$weights


ferman_scm <- ferman_for_analysis %>%
  get_cal_matches(
    covs = covs,
    treatment = "Z",
    caliper = c,
    metric = "maximum",   # "maximum", "euclidean", "manhattan"
    rad_method = "adaptive",
    dist_scaling = scaling,
    est_method = "scm",
    return = "sc_units",
    knn = 10,         # for ada
    num_bins = 5,     # for cem
    wider = F)
source("R/diagnostic_plots.R")
ferman_scm$id <-
  gsub("*_syn","",ferman_scm$id)
ferman_scm$id <- as.integer(ferman_scm$id)
# ferman_scm_adacalipers <-
#   attr(ferman_scm, "adacalipers")
# ferman_scm_adacalipers$id <-
#   as.character(ferman_scm_adacalipers$id)
# attr(ferman_scm, "adacalipers") <-
#   ferman_scm_adacalipers
love_plot(
  res = ferman_scm,
  covs = covs,
  B=NA) +
  ylim(c(-0.0065, 0.006))
ggsave(
  filename="writeup/figures/love-plot-ferman.png",
  width=8.76,
  height=5.33
  )


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
source("R/distance.R")
source("R/estimate.R")
ess_plot(feasible)


#####
## Get the balance table
#####
d <- feasible
metric <- "maximum"

sc_dists <-
  gen_dm(df = d %>%
           agg_sc_units(),
         covs = covs,
         treatment = "Z",
         scaling = scaling,
         metric = metric)%>%
  diag()

avg_dists <-
  gen_dm(df = d %>%
           agg_avg_units(),
         covs = covs,
         treatment = "Z",
         scaling = scaling,
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
         covs = covs,
         treatment = "Z",
         scaling = scaling,
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

## Number of used controls
n_t_SCM <- length(unique(
  (d %>% filter(Z==0, weights!=0))$id ))
n_t_avg <-length(unique(
  (d %>% filter(Z==0))$id ))
feasible_subclasses <- attr(ferman_scm, "feasible_subclasses")
n_feasible <- length(feasible_subclasses)
(n_t_1nn <- n_feasible)

n_c <- sum(ferman_for_analysis$Z==0)

######
# Histogram
######
matched_controls <-
  attr(ferman_scm, "scweights") %>%
  bind_rows() %>%
  group_by(subclass) %>%
  summarise(n_controls = n()-1) %>%
  filter(subclass %in% subclass_feasible)

sum(matched_controls$n_controls)
hist_matched_controls <-
  ggplot(matched_controls, aes(x = n_controls)) +
  geom_histogram(binwidth = 3, fill = "skyblue", color = "white") +
  labs(
    x = "Number of Matched Controls",
    y = "Frequency"
  ) +
  theme_minimal()

matched_controls_nonzero_weights <-
  attr(ferman_scm, "scweights") %>%
  bind_rows() %>%
  filter(weights > 0) %>%
  group_by(subclass) %>%
  summarise(n_controls = n()-1) %>%
  filter(subclass %in% subclass_feasible)
sum(matched_controls_nonzero_weights$n_controls)
hist_matched_controls_nonzero_weights<-
  ggplot(matched_controls_nonzero_weights, aes(x = n_controls)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "white") +
  labs(
    x = "Number of Matched Controls",
    y = "Frequency"
  ) +
  theme_minimal()
p <- grid.arrange(
  hist_matched_controls,
  hist_matched_controls_nonzero_weights,
  ncol=2
)
ggsave(
  "writeup/figures/hist-n-co.png",
  plot=p,
  width=8.1,
  height=5.3
)

#####
## SATT plot
######
source("R/diagnostic_plots.R")
ggd_att <- ferman_scm %>%
  left_join(attr(ferman_scm, "adacalipers"), by="id") %>%
  group_by(subclass) %>%
  summarize(adacal = last(adacal),
            tx = Y[2] - Y[1]) %>%
  arrange(adacal) %>%
  mutate(order = 1:n(),
         cum_avg = cumsum(tx) / order ) %>%
  slice((n_feasible):n())
plot_max_caliper_size <-
  ggd_att %>%
  ggplot(aes(x=order, y=adacal)) +
  geom_line(alpha=0.5) +
  geom_point(size=3) +
  theme_classic() +
  labs(y = "Maximum caliper size used",
       x = "Total number of treated units used") +
  expand_limits(color=1)


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
p
source("./R/bootstrap.R")

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

  # Calculate the statistic and save it to the vector
  se_AEs[i] <- get_se_AE(df_curr)
}

# Output the results
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

plot_all <- grid.arrange(
  plot_max_caliper_size,
  plot_SATT
)

ggsave(
  filename="writeup/figures/fsatt-ferman.png",
  plot=plot_all,
  width=5.3,
  height=7.7
)
