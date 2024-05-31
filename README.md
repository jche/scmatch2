
# Synthetic Caliper Matching (scmatch2)

This package implements the Synthetic Caliper Matching approach, which
is a blend of radius matching using distance metrics put on the
covariate distribution itself, and the synthetic control method. In
particular, it identifies sets of units local to each treated unit in
turn, and then makes a synthetic control for each treated unit using
those local units.

# Dependencies

This package has some tricky dependencies. In particular, it uses a DGP
(for simulation) from ACIC 2016:

    remotes::install_github("vdorie/aciccomp/2016")

You should not need it unless generating synthetic data.

# TO DO List

- Move the wrappers file to the scripts/simulation folder (it is
  utilities for simulation, not core package stuff)
- Keep boot_CSM (and the utility functions it uses), but move the other
  bootstrap stuff out of the main package.
- 

# Notes on scripts directory

In the scripts directory are several folders:

- analysis: This has the paper results. This makes plots and tables out
  of the simulation results.
- datagen: DGPs for demonstration. The simulation DGPs are in the
  sim_data.R file in the main package.
- demo: Illustrates use of the package.
- ferman-data-analysis: Empirical example in paper
- figs: Similar role as analysis, making figures out of existing results
- sims: Code that runs the various simulations in the paper

TODO: Reorg the directory structure a bit:

Thing we want:

- figs: Makes the illustration figures from the main paper, but nothing
  else. (No empirical or simulation results in this folder.)
- empirical: Does the empirical analysis (both Ferman data and LaLonde
  data analysis files)
- demo: Illustrates use of the package (only). Not anything included in
  the main paper.
- sims: The code to run the simulations
- sim_analysis: The code to analyze the simulation results
- datagen: Move this folder to old code, unless it is used by something
  else?

<!-- README.md is generated from README.Rmd. Please edit that file -->

# CSM

<!-- badges: start -->
<!-- badges: end -->

The goal of CSM is to …

## Installation

You can install the development version of CSM from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jche/scmatch2")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(CSM)
#> Loading required package: tidyverse
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.4     ✔ readr     2.1.5
#> ✔ forcats   1.0.0     ✔ stringr   1.5.1
#> ✔ ggplot2   3.5.0     ✔ tibble    3.2.1
#> ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
#> ✔ purrr     1.0.2     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
