
# functions to generate diagnostic plots




# SCM evaluation plot -----------------------------------------------------

dist_density_plot <- function(d, dist_scaling, metric) {
  # check distances bw each tx/sc pair
  sc_dists <- d %>% 
    agg_sc_units() %>% 
    gen_dm(scaling = dist_scaling,
           method = metric) %>% 
    diag()
  avg_dists <- d %>% 
    agg_avg_units() %>% 
    gen_dm(scaling = dist_scaling,
           method = metric) %>% 
    diag()
  
  nn_dists <- d %>% 
    group_by(Z,subclass) %>% 
    filter(dist == min(dist)) %>% 
    slice(1) %>% 
    mutate(weights = 1) %>% 
    ungroup() %>% 
    agg_avg_units() %>% 
    gen_dm(scaling = dist_scaling,
           method = metric) %>% 
    diag()
  
  tibble(SCM = sc_dists,
         Average = avg_dists,
         "1-NN" = nn_dists) %>% 
    pivot_longer(everything()) %>% 
    mutate(name = factor(name, levels=c("SCM", "Average", "1-NN"))) %>% 
    ggplot(aes(x=value, color=name)) +
    geom_density(linewidth=1) +
    theme_classic() +
    scale_color_manual(values = wesanderson::wes_palette("Zissou1", 5)[c(5,4,1)]) +
    labs(y = "Density",
         x = "Distance between treated unit and corresponding control",
         color = "Method used \nto generate \ncorresponding \ncontrol")
}

scm_vs_avg_plot <- function(d, dist_scaling, metric) {
  # check distances bw each tx/sc pair
  sc_dists <- d %>% 
    agg_sc_units() %>% 
    gen_dm(scaling = dist_scaling,
           method = metric) %>% 
    diag()
  avg_dists <- d %>% 
    agg_avg_units() %>% 
    gen_dm(scaling = dist_scaling,
           method = metric) %>% 
    diag()
  
  
  # identify places where average does better than sc
  #  (this should never happen)
  if (F) {
    ids <- feasible %>% 
      agg_sc_units() %>% 
      filter(!is.na(id)) %>% 
      pull(id)
    tibble(sc_dists = sc_dists,
           avg_dists = avg_dists) %>% 
      mutate(id = ids) %>% 
      filter(sc_dists > avg_dists) %>% 
      print(n=30)
    
    foo <- 153
    feasible %>% 
      filter(subclass == foo)
    feasible %>% 
      agg_sc_units() %>% 
      filter(subclass == foo)
    feasible %>% 
      agg_avg_units() %>% 
      filter(subclass == foo)
  }
  
  # TODO: simple plot comparing distance between:
  #  - tx unit and simple average control
  #  - tx unit and synthetic control
  
  tibble(sc_dists = sc_dists,
         avg_dists = avg_dists) %>% 
    ggplot(aes(x=sc_dists, y=avg_dists)) +
    geom_point() +
    geom_abline(lty="dotted") +
    theme_classic() +
    labs(y = "Distance between treated and average unit",
         x = "Distance between treated and SC unit")
}


# ESS plot ----------------------------------------------------------------

ess_plot <- function(d) {
  df_sc <- d %>% 
    filter(Z==F)
  
  df_avg <- d %>% 
    filter(Z==F) %>%
    group_by(subclass) %>%
    mutate(weights = 1/n()) %>%
    ungroup()
  
  df_1nn <- d %>% 
    filter(Z==F) %>% 
    group_by(subclass) %>% 
    filter(dist == min(dist)) %>% 
    slice(1) %>% 
    mutate(weights = 1) %>% 
    ungroup()
  
  df_tx <- d %>% 
    filter(Z==T)
  
  vars_df <- list(df_tx, df_sc, df_avg, df_1nn) %>% 
    map_dbl(~ .x %>% 
              agg_co_units() %>% 
              summarize(var = sum(weights^2)) %>% 
              pull(var)) %>% 
    tibble(method = c("tx", "sc", "avg", "1nn"),
           var = .)
  
  vars_df %>% 
    mutate(method = factor(method, levels=c("tx", "sc", "avg", "1nn"))) %>% 
    ggplot(aes(x=method, y=var)) +
    geom_col(aes(fill = method)) +
    scale_fill_manual(values = c("black",
                                 wesanderson::wes_palette("Zissou1", 5)[c(5,4,1)])) +
    scale_x_discrete(labels = c("Treated units",
                                "Control units \n(SCM weights)", 
                                "Control units \n(Average weights)", 
                                "Control units \n(1-NN weights)")) +
    theme_classic() +
    guides(fill="none") +
    labs(y = "Variance associated \nwith set of units",
         x = "")
}


# estimate-estimand tradeoff plot -----------------------------------------

satt_plot <- function(res, B=NA) {
  feasible_subclasses <- attr(res, "feasible_subclasses")
  n_feasible <- length(feasible_subclasses)
  
  ggd_maxcal <- res %>% 
    filter(!subclass %in% feasible_subclasses) %>% 
    filter(Z==1) %>% 
    left_join(attr(res, "adacalipers"), by="id") %>% 
    arrange(adacal) %>% 
    mutate(order = 1:n() + n_feasible)
  
  p1 <- ggd_maxcal %>% 
    ggplot(aes(x=order, y=adacal)) +
    geom_col(fill="gray") +
    geom_hline(yintercept=CALIPER, lty="dotted") +
    theme_classic() +
    labs(y = "Largest adaptive \ncaliper used",
         x = "Total number of treated units used")
  
  # ATT estimate vs. # co units added
  ggd_att <- res %>%
    left_join(attr(res, "adacalipers"), by="id") %>% 
    group_by(subclass) %>%
    summarize(adacal = last(adacal),
              tx = Y[2] - Y[1]) %>%
    arrange(adacal) %>%
    mutate(cum_avg = cumsum(tx) / (1:n()) ) %>% 
    slice((length(attr(res, "feasible_units"))+1):n()) %>% 
    mutate(order = 1:n() + n_feasible) 
  
  p2 <- ggd_att %>% 
    ggplot(aes(x=order, y=cum_avg)) +
    geom_point() +
    geom_step(linewidth=1) +
    theme_classic() +
    labs(y = "Cumulative ATT Estimate",
         x = "Total number of treated units used")
  
  if (!is.na(B)) {
    # bootstrap dists of FSATT and SATT
    boot_fsatt <- attr(res, "scweights") %>% 
      bind_rows() %>% 
      filter(subclass %in% feasible_subclasses) %>% 
      agg_co_units() %>% 
      boot_bayesian(B=B)
    
    boot_satt <- attr(res, "scweights") %>% 
      agg_co_units() %>% 
      boot_bayesian(B=B)
    
    # # using sd of bootstrap samples
    # boot_df <- ggd_att %>% 
    #   filter(order == min(order) | order == max(order)) %>% 
    #   mutate(sd = c(sd(boot_fsatt), sd(boot_satt)))
    # 
    # p2 <- p2 +
    #   geom_errorbar(data=boot_df,
    #                 aes(ymin=cum_avg-2*sd, ymax=cum_avg+2*sd),
    #                 width = 1)
    
    # using full distribution of bootstrap samples
    boot_df <- ggd_att %>% 
      filter(order == min(order) | order == max(order)) %>% 
      mutate(q025 = c(quantile(boot_fsatt, 0.025),
                      quantile(boot_satt, 0.025)),
             q975 = c(quantile(boot_fsatt, 0.975),
                      quantile(boot_satt, 0.975)))
    
    p2 <- p2 +
      geom_errorbar(data=boot_df,
                    aes(ymin=q025, ymax=q975),
                    width = 1)
  }
  
  p1/p2
}


