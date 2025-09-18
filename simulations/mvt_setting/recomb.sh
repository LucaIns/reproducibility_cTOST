#!/bin/sh
#SBATCH --job-name=cTOST_recomb
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-1:00:00
#SBATCH --partition=shared-cpu
#SBATCH --output simulations/sim_HPC_out/outfile/outfile_recomb_HPC.out

source ~/simulations/setting.sh

module load GCC/11.2.0 OpenMPI/4.1.1 R/4.2.0
srun Rscript --verbose $INFILE_recomb $sim_setting $n_simu > $OUTFILE_recomb