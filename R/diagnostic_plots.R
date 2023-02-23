
# functions to generate diagnostic plots

# compute full distance matrix for sc units (fast)
df_dm <- gen_dm(calada_scm,
                scaling = dm_scaling,
                method = METRIC)
sc_dists <- diag(df_dm)
calada_scm_dists <- calada_scm %>% 
  mutate(dist = rep(sc_dists, each=2))

# just look at some sc units
calada_scm_dists %>% 
  mutate(Y = round(Y, 3),
         dist = round(dist, 3)) %>% 
  arrange(dist, subclass, Z) %>%
  select(subclass, Z, X1, X2, X7, X8, Y, dist)

# check distances of all sc units
calada_scm_dists %>% 
  filter(Z==T) %>% 
  ggplot(aes(x=dist)) +
  geom_density(linewidth=1) +
  geom_vline(xintercept=CALIPER, lty="dashed") +
  theme_classic() +
  labs(y = "",
       x = "Scaled distance between tx and sc unit")