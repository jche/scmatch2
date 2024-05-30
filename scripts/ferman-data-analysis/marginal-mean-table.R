library(haven)
library(tidyverse)
library(xtable)

data_ferman <- read_dta(
  file = "./data/inputs/Final.dta"
)
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
         X3 = y2009)
# Make table: two columns: Z == 1 and Z ==0 represent
#  treated and control; three rows: X1, X2, X3
#   Each cell fill in average value of Y
table_result <- df_for_analysis %>%
  group_by(Z) %>%
  summarise(
    `Score 2007` = round(mean(X1, na.rm = TRUE), 6),
    `Score 2008` = round(mean(X2, na.rm = TRUE), 6),
    `Score 2009` = round(mean(X3, na.rm = TRUE), 6)
  ) %>%
  pivot_longer(cols = c(`Score 2007`,
                        `Score 2008`,
                        `Score 2009`),
               names_to = "Variable",
               values_to = "Average_Y") %>%
  pivot_wider(names_from = Z, values_from = Average_Y,
              names_prefix = "Z_") %>%
  rename(
    `Treated (Z == 1)` = Z_1,
    `Control (Z == 0)` = Z_0
  )

# Use R code to generate the latex code of table_result
latex_code <-
  xtable::xtable(table_result)
print(latex_code)

