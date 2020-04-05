#!/bin/zsh
#SBATCH -t 5
#Sbatch -n 1
#SBATCH -A clinical_analytics_lab

module load gcc/7.1.0
module load R/3.6.1

date

scriptfile=`Rscript -e 'system.file("scripts/run-mcmc.R", package="CovMitigation")'`

tid=$SLURM_ARRAY_TID
outfile="$1-$tid.rds"

echo "Rscript -e "

date