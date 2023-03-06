
# functions to generate diagnostic plots




# estimate-estimand tradeoff plot -----------------------------------------

satt_plot <- function(res, B=NA) {
  ggd_maxcal <- res %>% 
    filter(!subclass %in% feasible_subclasses) %>% 
    filter(Z==1) %>% 
    left_join(attr(res, "adacalipers"), by="id") %>% 
    arrange(adacal) %>% 
    mutate(order = 1:n())
  
  p1 <- ggd_maxcal %>% 
    ggplot(aes(x=order, y=adacal)) +
    geom_col(fill="gray") +
    geom_hline(yintercept=CALIPER, lty="dotted") +
    theme_classic() +
    labs(y = "Adaptive caliper",
         x = "Number of units added")
  
  # ATT estimate vs. # co units added
  ggd_att <- res %>%
    left_join(attr(res, "adacalipers"), by="id") %>% 
    group_by(subclass) %>%
    summarize(adacal = last(adacal),
              tx = Y[2] - Y[1]) %>%
    arrange(adacal) %>%
    mutate(cum_avg = cumsum(tx) / (1:n()) ) %>% 
    slice((length(feasible_units)+1):n()) %>% 
    mutate(order = 1:n()) 
  
  p2 <- ggd_att %>% 
    ggplot(aes(x=order, y=cum_avg)) +
    geom_point() +
    geom_step(linewidth=1) +
    theme_classic() +
    labs(y = "Cumulative ATT Estimate",
         x = "Number of units added")
  
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


