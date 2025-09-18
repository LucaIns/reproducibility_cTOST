#!/bin/sh

#SBATCH --job-name=cTOST_sims
#SBATCH --partition=shared-cpu
#SBATCH --time=0-08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output simulations/sim_HPC_out/outfile/outfile_%a.out

source ~/simulations/setting.sh

module load GCC/11.2.0 OpenMPI/4.1.1 R/4.2.0
srun Rscript --verbose $INFILE_sim $sim_setting $n_simu $n_array > $OUTFILE_sim