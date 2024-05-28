# Synthetic Caliper Matching (scmatch2)

This package implements the Synthetic Caliper Matching approach, which is a blend of radius matching using distance metrics put on the covariate distribution itself, and the synthetic control method.
In particular, it identifies sets of units local to each treated unit in turn, and then makes a synthetic control for each treated unit using those local units.





# Dependencies

This package has some tricky dependencies.
In particular, it uses a DGP (for simulation) from ACIC 2016:
```
remotes::install_github("vdorie/aciccomp/2016")
```
You should not need it unless generating synthetic data.
