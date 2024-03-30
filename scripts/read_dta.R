library(haven)

data_ferman <- read_dta(
  file = "./data/inputs/Final.dta"
)
colnames(data_ferman)
