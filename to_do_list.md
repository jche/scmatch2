

# TO DO List

- Move the wrappers file to the scripts/simulation folder (it is
                                                           utilities for simulation, not core package stuff)
- Keep boot_CSM (and the utility functions it uses), but move the other
bootstrap stuff out of the main package.
- Make a separate folder for the inference simulation and it should have
the annalysis file and the simulation file in it.

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
- sim_inference: The simulation and analysis code for the simulations
evaluating inference.
- datagen: Move this folder to old code, unless it is used by something

TODO:Improve functions in aggregation.R
- agg_sc_units should not take X, Y, Z
- agg_co_units should not set things to NA?

Improve function in get_cal_matches.R
- get_cal_matches takes one-dim covariates no scaling makes scaling variable not integer but $X value 1
