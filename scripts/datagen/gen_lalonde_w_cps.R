library(tidyverse)

lalonde_w_cps <- read_csv("data/inputs/lalonde.csv") %>%
  rename(Z = treat, Y = re78) %>%
  mutate(id = 1:n(),
         dataset = c(rep("NSWD", 185+260),
                     rep("CPS", 15992)),
         across(c(black, hispanic, married, nodegree), as.logical)) %>%
  select(id, Z, Y,
         black, hispanic, married, nodegree,
         age, education, re74, re75,
         dataset)

# rename variables for my matching methods
# xmeng: can rename at the matching's argument
# df <- lalonde %>%
#   rename(X1=black, X2=hispanic,
#          X3=married, X4=nodegree,
#          X5=age, X6=education,
#          X7=re74, X8=re75)

# option 1: "lalonde_only": use only experimental data
lalonde_df <- lalonde_w_cps %>% filter(dataset == "NSWD")
save(lalonde_df, file = "data/inputs/lalonde_only.RData")

# option 2: "lalonde_w_cps": use experimental treateds and all controls
lalonde_df <- lalonde_w_cps
save(lalonde_df, file = "data/inputs/lalonde_w_cps.RData")

# option 3: use experimental treateds and nonexperimental controls
lalonde_df <- lalonde_w_cps %>%
  filter((dataset == "NSWD" & Z==T) | (dataset=="CPS" & Z==F))
save(lalonde_df, file = "data/inputs/lalonde_w_cps_cos.RData")


# # EDA ---------------------------------------------------------------------
#
# df %>%
#   group_by(dataset, Z) %>%
#   summarize(across(starts_with("X"),
#                    mean))

if (F) {
  # NSWD is...
  #  - younger, less married
  #  - less educated, fewer degrees
  #  - much more black
  #  - much lower earnings
  df %>%
    group_by(dataset) %>%
    summarize(across(X1:X8, mean))

  # NSWD has...
  #  - less variance in age/educ/income
  df %>%
    group_by(dataset) %>%
    summarize(across(X1:X8, sd))

  # visualize covariate distributions
  df %>%
    group_by(dataset) %>%
    pivot_longer(X1:X8) %>%
    ggplot(aes(x=value, color=dataset)) +
    geom_density() +
    facet_wrap(~name, scales="free")
}

