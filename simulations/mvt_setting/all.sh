#!/bin/sh

source ~/simulations/setting.sh

ID=$(sbatch --parsable --array=1-$n_array ~/simulations/launch.sh)
ID=$(sbatch --parsable --dependency=afterany:${ID} ~/simulations/recomb.sh)
sbatch --dependency=afterok:${ID} ~/simulations/clean.sh