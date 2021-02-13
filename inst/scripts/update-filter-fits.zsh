#!/bin/zsh

#SBATCH -t 7200
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -A clinical_analytics_lab

module purge
module load gcc/7.1.0  openmpi/3.1.4
module load gsl/2.4
module load R/4.0.0

echo "start: " `date`

# Tasks should be 1-130
tid=$SLURM_ARRAY_TASK_ID

echo "Run command:"
echo Rscript -e "library(CovMitigation); doParallel::registerDoParallel(4); update_filter_models('./inputdata', $tid)"

Rscript -e "library(CovMitigation); doParallel::registerDoParallel(4); update_filter_models('./filter-updates.2020-09-06', $tid)"

echo "end: " `date`
