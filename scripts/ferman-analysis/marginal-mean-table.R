
source( here::here( "scripts/ferman-analysis/02-core-ferman-analysis.R" ) )
data_ferman = ferman_for_analysis

colnames(data_ferman)

df_for_analysis <-
  data_ferman %>%
  select(starts_with("y20"),
         Control,
         UF) %>%
  filter(if_all(starts_with("y20"), ~ !is.na(.x))) %>%
  filter(!is.na(Control)) %>%
  mutate(Y = y2010,
         Z = Control,
         X1 = y2007,
         X2 = y2008,
         X3 = y2009,
         X4 = UF == 35)


table_result <- df_for_analysis %>%
  group_by(Z) %>%
  summarise(
    `Score 2007_mean` = mean(X1, na.rm = TRUE),
    `Score 2007_sd`   = sd(X1, na.rm = TRUE),
    `Score 2007_n`    = sum(!is.na(X1)),

    `Score 2008_mean` = mean(X2, na.rm = TRUE),
    `Score 2008_sd`   = sd(X2, na.rm = TRUE),
    `Score 2008_n`    = sum(!is.na(X2)),

    `Score 2009_mean` = mean(X3, na.rm = TRUE),
    `Score 2009_sd`   = sd(X3, na.rm = TRUE),
    `Score 2009_n`    = sum(!is.na(X3)),

    `Pct Sao Paulo_mean` = mean(X4, na.rm = TRUE),
    `Pct Sao Paulo_sd`   = sd(X4, na.rm = TRUE),
    `Pct Sao Paulo_n`    = sum(!is.na(X4)),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -Z,
    names_to = c("Variable", ".value"),
    names_pattern = "^(.*)_(mean|sd|n)$"
  ) %>%
  mutate(
    mean = round(mean, 6),
    sd   = round(sd, 6)
  ) %>%
  pivot_wider(
    names_from = Z,
    values_from = c(mean, sd, n),
    names_glue = "{.value}_Z_{Z}"
  )

table_result <- table_result %>%
  mutate( std_diff =
            ( mean_Z_1 - mean_Z_0 ) /
            sqrt( ( sd_Z_1^2 + sd_Z_0^2 ) / 2 )
  )
table_result

if ( FALSE ) {
  table_result <- table_result %>%
  rename(
    `Treated mean (Z == 1)`  = mean_Z_1,
    `Control mean (Z == 0)`  = mean_Z_0,
    `Treated sd (Z == 1)`    = sd_Z_1,
    `Control sd (Z == 0)`    = sd_Z_0,
    `Treated n (Z == 1)`     = n_Z_1,
    `Control n (Z == 0)`     = n_Z_0
  ) %>%
  mutate(
    `Difference (treated - control)` =
      `Treated mean (Z == 1)` - `Control mean (Z == 0)`
  ) %>%
  select(
    Variable,
    `Treated n (Z == 1)`, `Treated mean (Z == 1)`, `Treated sd (Z == 1)`,
    `Control n (Z == 0)`, `Control mean (Z == 0)`, `Control sd (Z == 0)`,
    `Difference (treated - control)`
  )

}

table_result

# generate the latex code of table_result
latex_code <-
  xtable::xtable(table_result)
print(latex_code, include.rownames = FALSE)

