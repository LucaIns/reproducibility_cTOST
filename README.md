# reproducibility_cTOST

This repository contains the following files:

- `run_case_study.R`: to reproduce the results of the case study

- `run_simulations.R`: to reproduce all simulation results

These two scripts can be re-run as they are to obtain the tables and figures presented in the main text and Supplementary Material. 
By default, they re-use pre-computed intermediate results that may take long to run and/or require HPC.
Specifically, these files source the following folders:

- `case_study`: contains the code to reproduce the analysis of the case study presented in the paper, as well as the raw data

- `simulations`: contains the code to reproduce all simulation results presented in the paper (see instructions below)

- `aux_fun`: all functions used in the univariate and multivariate settings (please note that these functions are not supported anymore and kept only for reproducibility purposes, use the latest cTOST R package for analysis: https://github.com/stephaneguerrier/cTOST)

## Simulation Instructions

To re-run the simulations, we provide an example for `reproducibility_cTOST\simulations\mvt_setting` which extends to all of our code run on the HPC:

- if the code is run locally, just set `clusterUse=T` at the beginning of `sim_HPC_emp_power_func_multiv.R`

- if the code is run on the HPC, keep `clusterUse=F` at the beginning of `sim_HPC_emp_power_func_multiv.R` and then on the HPC:

	- set your path where the main `reproducibility_cTOST` folder (i.e. the GitHub repo transferred on the cluster) is located

	- run `sbatch reproducibility_cTOST/simulations/mvt_setting/all.sh`, this will in turn:

		1. source the settings in `setting.sh`

		2. source `launch.sh` to launch an array job based on `n_array` arrays to perform `n_simu` simulations (note: `n_array <= n_simu`) according to the `INFILE_sim` R script. This will store: `n_simu` simulation results in `reproducibility_cTOST/simulations/mvt_setting/sim_HPC_out/sim_setting_X/tmp`, as well as log files in `reproducibility_cTOST/simulations/mvt_setting/sim_HPC_out/sim_setting_$sim_setting/report/` and `reproducibility_cTOST/simulations/mvt_setting/sim_HPC_out/outfile/` (the latter are only useful in case of errors -- e.g., no output produced -- or to show R console outputs)
    
		3. as the step above is completed, it sources `recomb.sh` to launch a single job that iteratively opens and combines all file in `./tmp` through the `INFILE_recomb` R script. This will create an Rdata file in `reproducibility_cTOST/simulations/mvt_setting/sim_HPC_out/sim_setting_X` for each method used in the `INFILE_sim` simulation

		4. as the step above is completed, it sources `clean.sh` to clean all temporary and log files (remove this line if logs have to be checked)
  
- Finally, once all simulation results are saved, one can locally run `reproducibility_cTOST/simulations/mvt_setting/plots.R` to reproduce the plots and store them in a `tex` folder (which by default re-uses saved data rather than re-running all simulations)
