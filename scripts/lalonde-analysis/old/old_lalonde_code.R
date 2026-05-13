
# Define 'ess' function (must be here before the loop)
ess <- function( weights ) {
  sum( weights )^2 / sum( weights^2 )
}

# 1. Prepare the full data and Calculate Cumulative Stats GLOBALLY
# We must calculate order and cum_avg on the FULL set first,
# so that the "first" point in our plot represents the FSATT + 1 hard unit.
ggd_att_full <- result_table(lalonde_scm) %>%
  left_join(lalonde_scm$treatment_table, by = c("subclass")) %>%
  group_by(subclass) %>%
  summarize(
    adacal = last(na.omit(adacal)),
    tx = Y[Z==1] - mean(Y[Z==0])
  ) %>%
  ungroup() %>%
  arrange(adacal) %>% # Important: Sort by difficulty
  mutate(
    # This order represents the total number of treated units used
    order = 1:n(),
    # Cumulative ATT Estimate for the set 1:n
    cum_avg = cumsum(tx) / order
  )

# 2. Identify the "Base" set (Feasible units)
# These units are ALWAYS included in the calculation
base_feasible_subclasses <- ggd_att_full %>%
  filter(adacal <= CALIPER) %>%
  pull(subclass)

# 3. Create the plotting data (The "Hard" units)
# We filter for units > CALIPER, but keep the 'cum_avg' calculated from the full set
ggd_att_filtered <- ggd_att_full %>%
  filter(adacal > CALIPER)

# 4. Prepare for SE Loop
n_filtered <- nrow(ggd_att_filtered)
se_estimate_ATT_filtered <- numeric(n_filtered)

# Full unit table for pulling raw data
all_units_for_se <- result_table(lalonde_scm) %>%
  left_join(lalonde_scm$treatment_table, by = c("subclass","id")) %>%
  group_by(subclass) %>%
  mutate(adacal = ifelse(is.na(adacal), first(na.omit(adacal)), adacal)) %>%
  ungroup() %>%
  arrange(adacal)

# 5. Loop: Calculate SEs using (Base Set + Hard Units 1:i)
for (i in 1:n_filtered) {
  # i <- 1
  # A. Identify the specific "Hard" subclasses added at this step
  current_hard_subclasses <- ggd_att_filtered$subclass[1:i]

  # B. COMBINE Base Set + Current Hard Set
  # This ensures i=1 includes ALL feasible units plus the 1st hard unit
  combined_subclasses <- c(base_feasible_subclasses, current_hard_subclasses)

  # C. Filter the raw data to this combined set
  df_curr <- all_units_for_se %>%
    filter(subclass %in% combined_subclasses)

  # Check for valid variance groups (at least one control group with n>=2)
  has_valid_variance_group <- df_curr %>%
    filter(Z == 0) %>%
    group_by(subclass) %>%
    summarise(count = n()) %>%
    filter(count >= 2) %>%
    nrow() > 0

  if (!has_valid_variance_group) {
    se_estimate_ATT_filtered[i] <- NA
    next
  }

  # D. TOTAL VARIANCE ESTIMATION
  # We apply the V estimator to the cumulative dataset
  estimate_ATT_result <- estimate_ATT(
    matches = df_curr,
    outcome = "Y",
    treatment = "Z",
    variance_method = "pooled"
  )

  if (is.na(estimate_ATT_result$SE)) {
    se_estimate_ATT_filtered[i] <- NA
  } else {
    se_estimate_ATT_filtered[i] <- estimate_ATT_result$SE
  }
}

# 6. Plotting
plot_max_cal <- ggd_att_filtered %>%
  ggplot(aes(x=order, y=adacal)) +
  geom_line(alpha=0.5) +
  geom_point(size=1) +
  theme_classic() +
  labs(y = "Maximum caliper size used", x = "Total number of treated units used")

p_satt <- ggd_att_filtered %>%
  ggplot(aes(x=order, y=cum_avg)) +
  geom_line(alpha=0.5) +
  geom_point(size=1) +
  geom_hline(yintercept=0, lty="dotted") +
  geom_errorbar(
    aes(
      ymin = cum_avg - 1.96 * se_estimate_ATT_filtered,
      ymax = cum_avg + 1.96 * se_estimate_ATT_filtered
    ),
    width = 0.5, alpha = 0.2
  ) +
  labs(y = "Cumulative ATT Estimate", x = "Total number of treated units used") +
  theme_classic()

