#!/bin/zsh
#SBATCH -t 120
#SBATCH -n 2
#SBATCH -c 1
#SBATCH -A clinical_analytics_lab

module load intel/18.0
module load intelmpi/18.0
module load R/3.6.3

echo "start:  " `date`

program="./run-filter-fit-loc.R"

# Tasks should be 1-133
tid=$SLURM_ARRAY_TASK_ID

echo "Run command:"
echo "source('$program'); run_filter_fit_loc($tid)"

Rscript -e "source('$program'); run_filter_fit_loc($tid)"

echo "end:  " `date`
