
<!-- badges: start -->
<!-- badges: end -->

# Caliper Synthetic Matching (scmatch2)

This package implements the Caliper Synthetic Matching approach, which
is a blend of radius matching using distance metrics put on the
covariate distribution itself, and the synthetic control method. In
particular, it identifies sets of units local to each treated unit in
turn, and then makes a synthetic control for each treated unit using
those local units.

Details can be found on the paper: 

* **Original CSM Method:** [Che et. al. (2024), Caliper Synthetic Matching](https://arxiv.org/abs/2411.05246)
* **Variance Estimation (Paper 2):** [Preprint: Inference in Matching (November 2025)](https://mengeks.github.io/files/Inference_in_Matching-10Nov2025.pdf)

The GitHub repo also has replication materials for the associated paper
in the `scripts` directory; see below for further discussion of this
added-on code. The installed package will ignore `scripts`.

## Installation

You can install the development version of CSM from
[GitHub](https://github.com/) with:

    # install.packages("devtools")  # (if needed)
    devtools::install_github("jche/scmatch2")

# Quick demo of package

We generate a small toy dataset to illustrate the main methods of
interest:

``` r
set.seed( 4044440 )
dat <- gen_one_toy(nt = 5)
dat
#> # A tibble: 503 × 8
#>       id    X1    X2 Z       noise    Y0    Y1     Y
#>    <int> <dbl> <dbl> <lgl>   <dbl> <dbl> <dbl> <dbl>
#>  1     1 0.189 0.219 TRUE   0.365   5.41  6.64  6.64
#>  2     2 0.220 0.176 TRUE   0.163   5.19  6.38  6.38
#>  3     3 0.776 0.641 TRUE  -0.558   4.50  8.76  8.76
#>  4     4 0.660 0.644 TRUE  -0.324   4.91  8.82  8.82
#>  5     5 0.778 0.129 FALSE  0.593   3.72  6.44  3.72
#>  6     6 0.702 0.216 FALSE -0.105   3.84  6.59  3.84
#>  7     7 0.859 0.306 FALSE -0.207   3.40  6.90  3.40
#>  8     8 0.820 0.302 FALSE -0.455   3.33  6.70  3.33
#>  9     9 0.764 0.248 FALSE -0.594   3.21  6.24  3.21
#> 10    10 0.828 0.183 FALSE  0.0198  3.17  6.20  3.17
#> # ℹ 493 more rows
ggplot( dat, aes( X1, X2, color = Z ) ) + geom_point() +
  coord_fixed()
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="75%" style="display: block; margin: auto;" />

To calculate matches, call `get_cal_matches()`–it will match, make
synthetic controls for each unit, and give you a final dataset back,
stored as an `csm_matches` object:

``` r
mtch <- get_cal_matches( dat, 
                        metric = "maximum",
                        scaling = c( 1/0.2, 1/0.2),
                        caliper = 1, 
                        rad_method = "adaptive", 
                        est_method = "scm" ) 
mtch
#> csm_matches: matching w/ maximum distance on 2, 3 
#> 4 Treated units matched to control units (0 above set caliper) 
#> Adaptive calipers: 1, 1, 1, 1 
#>  Target caliper = 1 
#>  Max distance ranges 0.945 - 0.988 
#> scaling: 5, 5
```

There are a variety of things you can pull from the result object. You
can get a list of statistics on the treated units:

``` r
mtch$treatment_table
#> # A tibble: 4 × 8
#>   id    subclass    nc   ess max_dist adacal feasible matched
#>   <chr> <chr>    <dbl> <dbl>    <dbl>  <dbl>    <dbl>   <dbl>
#> 1 1     1           24  2.69    0.945      1        1       1
#> 2 2     2           22  2.08    0.970      1        1       1
#> 3 3     3           37  2.80    0.987      1        1       1
#> 4 4     4           44  2.18    0.988      1        1       1
```

You can see all the units used, grouped by subclass:

``` r
full_unit_table(mtch, nonzero_weight_only = TRUE ) 
#> # A tibble: 16 × 12
#>    id        X1     X2 Z       noise    Y0    Y1     Y  dist subclass unit 
#>    <chr>  <dbl>  <dbl> <lgl>   <dbl> <dbl> <dbl> <dbl> <dbl> <chr>    <chr>
#>  1 1     0.189  0.219  TRUE   0.365   5.41  6.64  6.64 0     1        tx1  
#>  2 373   0.0688 0.0655 FALSE -0.148   4.63  5.04  4.63 0.768 1        c4   
#>  3 374   0.0183 0.380  FALSE -0.493   3.79  4.99  3.79 0.856 1        c5   
#>  4 446   0.363  0.287  FALSE -0.0208  5.16  7.11  5.16 0.866 1        c18  
#>  5 2     0.220  0.176  TRUE   0.163   5.19  6.38  6.38 0     2        tx1  
#>  6 353   0.414  0.303  FALSE  0.474   5.64  7.79  5.64 0.970 2        c3   
#>  7 373   0.0688 0.0655 FALSE -0.148   4.63  5.04  4.63 0.755 2        c4   
#>  8 502   0.0409 0.315  FALSE -0.658   3.90  4.97  3.90 0.894 2        c21  
#>  9 3     0.776  0.641  TRUE  -0.558   4.50  8.76  8.76 0     3        tx1  
#> 10 16    0.716  0.510  FALSE  0.599   5.60  9.27  5.60 0.659 3        c2   
#> 11 357   0.733  0.836  FALSE  0.403   5.41 10.1   5.41 0.974 3        c15  
#> 12 467   0.971  0.491  FALSE  0.0523  3.91  8.30  3.91 0.976 3        c30  
#> 13 4     0.660  0.644  TRUE  -0.324   4.91  8.82  8.82 0     4        tx1  
#> 14 21    0.672  0.525  FALSE -0.191   4.95  8.54  4.95 0.593 4        c3   
#> 15 385   0.807  0.803  FALSE -0.207   4.83  9.66  4.83 0.794 4        c21  
#> 16 498   0.563  0.842  FALSE -0.957   3.75  7.96  3.75 0.988 4        c43  
#> # ℹ 1 more variable: weights <dbl>
```

You can get the final dataset (here filtering to only matches within the
initial caliper):

``` r
result_table( mtch, feasible_only = TRUE )
#> # A tibble: 8 × 9
#>   id    subclass Z        X1    X2    Y0    Y1     Y weights
#>   <chr> <chr>    <lgl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>
#> 1 1_syn 1        FALSE 0.189 0.219  4.71  5.94  4.71       1
#> 2 1     1        TRUE  0.189 0.219  5.41  6.64  6.64       1
#> 3 2_syn 2        FALSE 0.220 0.176  5.06  6.24  5.06       1
#> 4 2     2        TRUE  0.220 0.176  5.19  6.38  6.38       1
#> 5 3_syn 3        FALSE 0.776 0.641  5.17  9.42  5.17       1
#> 6 3     3        TRUE  0.776 0.641  4.50  8.76  8.76       1
#> 7 4_syn 4        FALSE 0.660 0.644  4.61  8.52  4.61       1
#> 8 4     4        TRUE  0.660 0.644  4.91  8.82  8.82       1
```

You can also get the final generated result as a data.frame by casting
the result into a dataframe:

``` r
head( as.data.frame( mtch ), n = 10 )
#> # A tibble: 8 × 9
#>   id    subclass Z        X1    X2    Y0    Y1     Y weights
#>   <chr> <chr>    <lgl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>
#> 1 1_syn 1        FALSE 0.189 0.219  4.71  5.94  4.71       1
#> 2 1     1        TRUE  0.189 0.219  5.41  6.64  6.64       1
#> 3 2_syn 2        FALSE 0.220 0.176  5.06  6.24  5.06       1
#> 4 2     2        TRUE  0.220 0.176  5.19  6.38  6.38       1
#> 5 3_syn 3        FALSE 0.776 0.641  5.17  9.42  5.17       1
#> 6 3     3        TRUE  0.776 0.641  4.50  8.76  8.76       1
#> 7 4_syn 4        FALSE 0.660 0.644  4.61  8.52  4.61       1
#> 8 4     4        TRUE  0.660 0.644  4.91  8.82  8.82       1
```

You can estimate impacts on the matched dataset using whatever tools you
want, or use the options discussed in our paper as so:

``` r
get_ATT_estimate( mtch )
#> # A tibble: 1 × 6
#>     ATT    SE sigma_hat   N_T N_C_tilde     t
#>   <dbl> <dbl>     <dbl> <int>     <dbl> <dbl>
#> 1  2.76 0.335     0.543     4      7.71  8.25
```

# Dependencies

This package has some tricky dependencies. In particular, it
(optionally) uses a DGP (for simulation) from ACIC 2016:

    remotes::install_github("vdorie/aciccomp/2016")

You should not need this package unless you are generating synthetic
data.
