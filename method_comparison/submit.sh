#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=correlation
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=6000 
#SBATCH --time=24:00:00


module load GSL
module load MPFR
module load R

n_topic=50

Rscript kmeans.R $n_topic 
#Rscript kmeans_corr.R 
