#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=tune
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=6000 
#SBATCH --time=12:00:00


module load GSL
module load MPFR
module load R

n_start=125
n_end=140
step=5

Rscript tuning_nonpara.r $n_start $n_end $step
