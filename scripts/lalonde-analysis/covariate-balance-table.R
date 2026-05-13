# scripts/lalonde-analysis/covariate-balance-table.R
# Table of covariate means and SDs by treatment/control group

library(tidyverse)
library(here)

lalonde_clean <- readRDS(here::here("scripts/lalonde-analysis/data/lalonde_for_analysis.rds"))

cov_labels <- c(
  X1 = "Black",
  X2 = "Hispanic",
  X3 = "Married",
  X4 = "No degree",
  X5 = "Age",
  X6 = "Education",
  X7 = "Earnings 1974",
  X8 = "Earnings 1975"
)

cov_table <- lalonde_clean |>
  pivot_longer(starts_with("X"), names_to = "covariate", values_to = "value") |>
  group_by(Z, covariate) |>
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value,   na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from  = Z,
    values_from = c(mean, sd),
    names_glue  = "{ifelse(Z == 1, 'Treated', 'Control')}_{.value}"
  ) |>
  mutate(covariate = cov_labels[covariate]) |>
  rename(Covariate = covariate) |>
  select(Covariate,
         `Treated_mean`, `Treated_sd`,
         `Control_mean`, `Control_sd`) |>
  mutate(SMD = (Treated_mean - Control_mean) / Treated_sd)

cov_table
