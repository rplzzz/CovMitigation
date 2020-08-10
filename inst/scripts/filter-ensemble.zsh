#!/bin/zsh
#SBATCH -t 2500
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -A clinical_analytics_lab

module load intel/18.0
module load intelmpi/18.0
module load R/3.6.3

echo "start:  " `date`

program="./gen-filter-ensemble.R"


## The valid task ids are 1-133 for the 133 localities for which we have data.
tid=$SLURM_ARRAY_TASK_ID

echo "Run command:"
echo "source('$program'); gen_ensemble($tid)"

Rscript -e "source('$program'); gen_ensemble($tid)"

echo "end:  " `date`
