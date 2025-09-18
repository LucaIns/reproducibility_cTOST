export sim_setting=1  # setting number (input: from 1 to 4)
export n_simu=50000   # number of simulations
export n_array=5000   # number of array for HPC

export INFILE_sim=simulations/sim_HPC_emp_power_func_multiv.R 
export OUTFILE_sim=simulations/sim_HPC_out/sim_setting_${sim_setting}/report/report_${SLURM_ARRAY_TASK_ID}.Rout

export INFILE_recomb=simulations/out_HPC_emp_power_func_multiv.R 
export OUTFILE_recomb=simulations/sim_HPC_out/sim_setting_${sim_setting}/report_recomb_HPC.Rout
