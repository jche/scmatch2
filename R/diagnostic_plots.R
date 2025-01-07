library(latex2exp)


# functions to generate diagnostic plots
create_toy_df_plot <- function(toy_df) {
  toy_df_plot <- toy_df %>%
    mutate(Z = factor(Z, levels = c(FALSE, TRUE), labels = c("Control", "Treatment"))) %>%
    ggplot(aes(x = X1, y = X2, col=as.factor(Z))) +
    geom_point() +
    # scale_color_continuous(low = "blue", high = "orange") +
    theme_classic() +
    labs(x = TeX("$X_1$"),
         y = TeX("$X_2$"),
         color = "Group")
  return(toy_df_plot)
}




# SCM evaluation plot -----------------------------------------------------

dist_density_plot <- function(d, scaling, metric) {
  # check distances bw each tx/sc pair
  sc_dists <- d %>%
    agg_sc_units() %>%
    gen_dm(scaling = scaling,
           metric = metric) %>%
    diag()
  avg_dists <- d %>%
    agg_avg_units() %>%
    gen_dm(scaling = scaling,
           metric = metric) %>%
    diag()

  nn_dists <- d %>%
    group_by(Z,subclass) %>%
    filter(dist == min(dist)) %>%
    slice(1) %>%
    mutate(weights = 1) %>%
    ungroup() %>%
    agg_avg_units() %>%
    gen_dm(scaling = scaling,
           metric = metric) %>%
    diag()

  res_list <- list(sc = sc_dists, avg = avg_dists, nn = nn_dists)
  print(tibble(
    method = c("SCM", "Average", "1-NN"),
    mean = map_dbl(res_list, ~round(mean(.), 3)),
    median = map_dbl(res_list, ~round(median(.), 3))
  ))

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
         x = "Distance between treated unit \nand corresponding control",
         color = "Method")
}

scm_vs_avg_plot <- function(d, scaling, metric) {
  # check distances bw each tx/sc pair
  sc_dists <- d %>%
    agg_sc_units() %>%
    gen_dm(scaling = scaling,
           metric = metric) %>%
    diag()
  avg_dists <- d %>%
    agg_avg_units() %>%
    gen_dm(scaling = scaling,
           metric = metric) %>%
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

#' Effective sample size (ESS) plot
#'
#' @export
ess_plot <- function(d) {
  df_sc <- d %>%
    filter(Z==FALSE)

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
              summarize(ess = sum(weights)^2 / sum(weights^2)) %>%
              pull(ess)) %>%
    tibble(method = c("tx", "sc", "avg", "1nn"),
           ess = .)

  print(vars_df)

  vars_df %>%
    mutate(method = factor(method, levels=c("1nn", "avg", "sc", "tx"))) %>%
    ggplot(aes(x=method, y=ess)) +
    geom_col(aes(fill = method)) +
    scale_fill_manual(values = c("black",
                                 wesanderson::wes_palette("Zissou1", 5)[c(5,4,1)])) +
    scale_x_discrete(labels = c("Control units \n(1-NN weights)",
                                "Control units \n(Average weights)",
                                "Control units \n(SCM weights)",
                                "Treated units")) +
    theme_classic() +
    guides(fill="none") +
    labs(y = "ESS of set of units",
         x = "") +
    coord_flip()
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

satt_plot2 <- function(res, B=NA) {
  feasible_subclasses <- attr(res, "feasible_subclasses")
  n_feasible <- length(feasible_subclasses)

  # ATT estimate vs. # co units added
  ggd_att <- res %>%
    left_join(attr(res, "adacalipers"), by="id") %>%
    group_by(subclass) %>%
    summarize(adacal = last(adacal),
              tx = Y[2] - Y[1]) %>%
    arrange(adacal) %>%
    mutate(order = 1:n(),
           cum_avg = cumsum(tx) / order ) %>%
    slice((n_feasible+1):n())

  foo <- ggd_att %>%
    group_by(adacal) %>%
    summarize(cum_avg = last(cum_avg),
              n = n()) %>%
    mutate(shape = ifelse(n==1, 20, 19))
  p <- ggd_att %>%
    ggplot(aes(x=adacal, y=cum_avg)) +
    geom_point(data = foo,
               aes(shape = shape),
               size=2) +
    geom_step() +
    theme_classic() +
    labs(y = "Cumulative ATT Estimate",
         x = "Maximum caliper size used") +
    scale_shape_identity()

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

    p <- p +
      geom_errorbar(data=boot_df,
                    aes(ymin=q025, ymax=q975),
                    width = 0.05)
  }

  p
}

satt_plot3 <- function(res, B=NA) {
  feasible_subclasses <- attr(res, "feasible_subclasses")
  n_feasible <- length(feasible_subclasses)

  # ATT estimate vs. # co units added
  ggd_att <- res %>%
    left_join(attr(res, "adacalipers"), by="id") %>%
    group_by(subclass) %>%
    summarize(adacal = last(adacal),
              tx = Y[2] - Y[1]) %>%
    arrange(adacal) %>%
    mutate(order = 1:n(),
           cum_avg = cumsum(tx) / order ) %>%
    slice((n_feasible):n())

  p <- ggd_att %>%
    ggplot(aes(x=order, y=cum_avg)) +
    geom_line(alpha=0.5) +
    geom_point(aes(color=adacal),
               size=3) +
    theme_classic() +
    labs(y = "Cumulative ATT Estimate",
         x = "Total number of treated units used",
         color = "Maximum caliper size used") +
    scale_color_continuous(low="blue", high="orange") +
    expand_limits(color=1)

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

    # using sd of bootstrap samples
    boot_df <- ggd_att %>%
      filter(order == min(order) | order == max(order)) %>%
      mutate(sd = c(sd(boot_fsatt), sd(boot_satt)))

    print(paste0("FSATT: (",
                 round(boot_df$cum_avg[1]-1.96*boot_df$sd[1],3),
                 ", ",
                 round(boot_df$cum_avg[1]+1.96*boot_df$sd[1], 3), ")"))
    print(paste0("SATT: (",
                 round(boot_df$cum_avg[2]-1.96*boot_df$sd[2], 3),
                 ", ",
                 round(boot_df$cum_avg[2]+1.96*boot_df$sd[2], 3), ")"))

    p <- p +
      geom_errorbar(data=boot_df,
                    aes(ymin=cum_avg-1.96*sd, ymax=cum_avg+1.96*sd),
                    width = 0.5,
                    linewidth=1) +
      geom_point(aes(color=adacal),
                 size=3)

    # # using full distribution of bootstrap samples
    # boot_df <- ggd_att %>%
    #   filter(order == min(order) | order == max(order)) %>%
    #   mutate(q025 = c(quantile(boot_fsatt, 0.025),
    #                   quantile(boot_satt, 0.025)),
    #          q975 = c(quantile(boot_fsatt, 0.975),
    #                   quantile(boot_satt, 0.975)))
    #
    # p <- p +
    #   geom_errorbar(data=boot_df,
    #                 aes(ymin=q025, ymax=q975),
    #                 width = 0.5) +
    #   geom_point(aes(color=adacal),
    #              size=3)
  }

  p
}

# fixed, with wild bootstrap
satt_plot4 <- function(res, B=NA) {
  feasible_subclasses <-
    attr(res, "feasible_subclasses")
  n_feasible <- length(feasible_subclasses)

  # ATT estimate vs. # co units added
  ggd_att <- res %>%
    left_join(attr(res, "adacalipers"), by="id") %>%
    group_by(subclass) %>%
    summarize(adacal = last(adacal),
              tx = Y[2] - Y[1]) %>%
    arrange(adacal) %>%
    mutate(order = 1:n(),
           cum_avg = cumsum(tx) / order ) %>%
    slice((n_feasible):n())

  p <- ggd_att %>%
    ggplot(aes(x=order, y=cum_avg)) +
    geom_line(alpha=0.5) +
    geom_point(
      # aes(color=adacal),
               size=3) +
    theme_classic() +
    labs(y = "Cumulative ATT Estimate",
         x = "Total number of treated units used",
         color = "Maximum caliper size used") +
    # scale_color_continuous(low="blue", high="orange") +
    expand_limits(color=1)

  if (!is.na(B)) {
    # bootstrap dists of FSATT and SATT
    boot_fsatt <- res %>%
      filter(subclass %in% feasible_subclasses) %>%
      wild_bootstrap(B=B)

    boot_satt <- attr(res, "scweights") %>%
      agg_sc_units() %>%
      wild_bootstrap(B=B)

    # using sd of bootstrap samples
    boot_df <- ggd_att %>%
      filter(order == min(order) | order == max(order)) %>%
      mutate(sd = c(sd(boot_fsatt), sd(boot_satt)))

    print(paste0("FSATT: (",
                 round(boot_df$cum_avg[1]-1.96*boot_df$sd[1],3),
                 ", ",
                 round(boot_df$cum_avg[1]+1.96*boot_df$sd[1], 3), ")"))
    print(paste0("SATT: (",
                 round(boot_df$cum_avg[2]-1.96*boot_df$sd[2], 3),
                 ", ",
                 round(boot_df$cum_avg[2]+1.96*boot_df$sd[2], 3), ")"))

    p <- p +
      geom_errorbar(data=boot_df,
                    aes(ymin=cum_avg-1.96*sd, ymax=cum_avg+1.96*sd),
                    width = 0.5,
                    linewidth=1) +
      geom_point(aes(color=adacal),
                 size=3)
  }

  p
}


# love plot ---------------------------------------------------------------

get_diff_scm_co_and_tx <- function(res,covs){
  stopifnot( is.csm_matches(res) )

  ada = res$treatment_table %>%
    select(id, adacal) %>%
    mutate( id = as.character(id) )
  df_diff_scm_co_and_tx <- result_table(res) %>%
    left_join(ada,
              by="id") %>%
    group_by(subclass) %>%
    summarize(adacal = last(adacal),
              across(all_of(covs),
                     ~.[2] - .[1]))
  return(df_diff_scm_co_and_tx)
}

create_love_plot_df <- function(res, covs){
  feasible_subclasses <- feasible_units(res)
  n_feasible <- nrow(feasible_subclasses)

  df_step_1<-
    get_diff_scm_co_and_tx(res,covs)

  df_step_2<-df_step_1 %>%
    arrange(adacal) %>%
    mutate(order = 1:n(),
           across(all_of(covs),
                  ~cumsum(.) / order))

  df_step_3 <- df_step_2 %>%
    slice((n_feasible):n()) %>%
    pivot_longer(all_of(covs))
  love_plot_df <- df_step_3
  return(love_plot_df)
}


#' Love plot of covariate balance
#'
#' Make a ggplot love plot of covariate balance for each covariate
#' passed.
#'
#' @export
love_plot <- function(res, covs, covs_names = NULL, B=NA) {
  love_steps <- create_love_plot_df(res, covs)
  if ( !is.null( covs_names ) ) {
    love_steps$name <- covs_names[ match( love_steps$name, covs ) ]
  }
  p <- love_steps %>%
    ggplot(aes(x=order, color=name)) +
    geom_point(data=. %>%
                 slice(c(1:length(covs),
                      (n()-length(covs)):n())),
               aes(y=value), size=2) +
    geom_step(aes(y=value, group=name),
              linewidth=1.1) +
    expand_limits(y = 0) +
    geom_hline(yintercept=0,
               lty="dotted") +
    facet_wrap(~name,
               scales="free_y") +
    labs(y = "\n Covariate balance (tx-co)",
         x = "Total number of treated units used",
         color = "Covariate") +
    theme_classic()

  if (!is.na(B)) {
    # bootstrap dists of FSATT and SATT
    boot_fsatt <- attr(res, "scweights") %>%
      bind_rows() %>%
      filter(subclass %in% feasible_subclasses) %>%
      agg_co_units() %>%
      boot_bayesian_covs(covs=covs, B=B) %>%
      pivot_longer(everything()) %>%
      group_by(name) %>%
      summarize(q025 = quantile(value, 0.025),
                q975 = quantile(value, 0.975)) %>%
      mutate(order = n_feasible+1)

    boot_satt <- attr(res, "scweights") %>%
      agg_co_units() %>%
      boot_bayesian_covs(covs=covs, B=B) %>%
      pivot_longer(everything()) %>%
      group_by(name) %>%
      summarize(q025 = quantile(value, 0.025),
                q975 = quantile(value, 0.975)) %>%
      mutate(order = max(love_steps$order))

    p <- p +
      geom_errorbar(data=bind_rows(boot_fsatt, boot_satt),
                    aes(ymin=q025, ymax=q975),
                    width = 1)
  }

  return(p)

}


love_plot2 <- function(res, covs, B=NA) {
  feasible_subclasses <- attr(res, "feasible_subclasses")
  n_feasible <- length(feasible_subclasses)

  adacal_key <- attr(res, "adacalipers") %>%
    rename(subclass = id) %>%
    arrange(adacal) %>%
    mutate(order = 1:n())

  love_steps <- res %>%
    left_join(adacal_key, by="subclass") %>%
    group_by(Z) %>%
    arrange(order) %>%
    mutate(across(all_of(covs), ~cumsum(.) / order)) %>%
    slice((n_feasible):n()) %>%
    pivot_longer(all_of(covs))

  love_steps %>%
    mutate(Z = as.factor(Z)) %>%
    ggplot(aes(x=order, y=value, color=Z)) +
    # geom_line(aes(group=order), alpha=0.3, color="black") +
    geom_point(size=1) +
    geom_point(data=. %>%
                 filter(order == max(order) |
                          order == min(order)),
               size=2) +
    geom_line(aes(group=interaction(name, Z)),
              alpha=0.3) +
    facet_wrap(~name, scales="free_y") +
    labs(y = "Marginal mean",
         x = "Total number of treated units used",
         color = "Z") +
    theme_classic()
}

