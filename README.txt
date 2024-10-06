# Decomposing dynamical subprocesses for compositional generalization

This is the code supplement for the statistical models/GLMs implemented in
R, simulations and fitting of successor feature models implemented in MATLAB

## Instructions

To run the code, you'll need to install RStudio (specific libraries: rethinking - uses RStan) and MATLAB (version 2023b).

Follow these steps to reproduce the results:
1.) Statistics/GLMs: use stats.Rmd in RStudio - preprocessed data (long data format) is loaded from subfolder "data". 
You will have to set paths to this folder in the script, depending on where it is stored on your computer.

2.) Model simulations (prior): use sim_prior_main.m in MATLAB to run prior simulations with successor feature/compound models (functions "*_model_all.m"). 
Scripts use task schedules from subfolder "schedules", schedule numbers from "sched_num_all.mat" and load parameter values from space_alpha.mat, space_omega.mat, space_tau.mat and 
Scripts use function "plotting simulations" from subfolder "functions" 
Scripts store results/plots in "prior_simulations". 

3.) Model fitting (MLE): use fit_models_main.m in MATLAB. Script calls "*_model_all.m" functions implementing the models and task schedules from 
subfolder "schedules". Uses helper functions from subfolder "functions", stores fitted parameters in subfolder "fitted_parameters" 
(also contains fitted parameters from Structural Inference model - fits_all_llmax.csv). 

Takes input arguments: 
-subj (subject number, 1:241)
-iter (which restart of fitting process, 1:20)
-model_type (models 1-5) 
-multi_lr (1 or 2, single or separate learning rate parameters for prior and transfer learning)
-multi_tau (1 or 2, single or separate softmax stochasticity parameters for experience and inference probes) 
-prior_sigma (needed if using MAP fitting, set to arbitrary value, e.g. 0.1)
-MAP_estimation (needed if using MAP fitting, 0 otherwise - for MLE)
-oracle (fit SF transfer model with "oracle", mixing in the true TM instead of prior learning SR)
-run_cluster (0 if running on local machine - will take a long time, or parallelized on cluster)

4.) Model comparison: use model_comparison.m in MATLAB. Loads fitted parameters, likelihood values and subject IDs/condition labels stored in 
subfolder "fitted_parameters". Plots BIC values and model frequencies.

5.) Model simulations (posterior): use sim_posterior_main.m in MATLAB to posterior simulations of successor feature models with fitted parameter estimates 
(loaded from subfolder "fitted_parameters"). Calls model functions ("*_model_all.m"). 
Stores results in subfolder "posterior_simulations". Needs parallel pool, iterations (1000) may take some time, depending on how many parallel workers you can call.

6.) Plot posterior simulations using "sim_plots_posterior.m" in MATLAB: Uses "all_sim_plots.m" and functions from subfolder "functions"
