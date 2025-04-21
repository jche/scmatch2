# install.packages("devtools")  # (if needed)
# devtools::install_github("jche/scmatch2")

library(CSM)

set.seed( 4044440 )
dat <- gen_one_toy(nt = 50)
dat

mtch <- get_cal_matches( dat,
                         metric = "maximum",
                         scaling = c( 1/0.2, 1/0.2),
                         caliper = 1,
                         rad_method = "adaptive",
                         est_method = "scm" )
mtch

mtch$treatment_table # a table of statistics on the treated units

full_unit_table(mtch, nonzero_weight_only = TRUE )
# table of all units used, grouped by subclass (cluster)
# of treated and match controls.
# There is also "dist": the distance between the treated and the control
# "weights": the weights assigned to the control units

result_table( mtch, feasible_only = TRUE )

# Effect estimate
get_ATT_estimate( mtch )
