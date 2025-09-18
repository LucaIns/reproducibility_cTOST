#!/bin/sh

#SBATCH --job-name=cTOST_clean
#SBATCH --partition=debug-cpu
#SBATCH --time=0-00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output simulations/sim_HPC_out/outfile/clean.out

source ~/simulations/setting.sh

rm simulations/sim_HPC_out/sim_setting_$sim_setting/tmp/*
rm simulations/sim_HPC_out/sim_setting_$sim_setting/report/*
rm simulations/sim_HPC_out/outfile/*