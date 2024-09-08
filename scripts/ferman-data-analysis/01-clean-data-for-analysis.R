library(tidyverse)
data_ferman <- read_dta(
  file = here::here( "data/inputs/Final.dta" )
)

colnames(data_ferman)
hist(data_ferman$y2008 - data_ferman$OBJETIVA2008)
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
         is_sao_paolo = as.numeric(UF == 35))
# UF == 35 means Sao Paolo
# UF == 33 means Rio

saveRDS(ferman_for_analysis,
        file = here::here( "data/inputs/ferman_for_analysis.rds" ))

# EDA
if (FALSE){
  proportion_by_UF <- ferman_for_analysis %>%
    group_by(Control, is_sao_paolo) %>%
    summarise(n())
}

