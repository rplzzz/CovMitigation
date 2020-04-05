#!/bin/zsh
#SBATCH -t 1000
#Sbatch -n 1
#SBATCH -A clinical_analytics_lab

module load gcc/7.1.0
module load R/3.6.1

echo "start:  " `date`

program=`Rscript -e 'system.file("scripts/run-mcmc.R", package="CovMitigation")'`

tid=$SLURM_ARRAY_TID
nsamp=100000
outfile="$1-$tid.rds"

echo "Run command:"
echo "source('$program'); run_mcmc($tid, $nsamp, '$outfile', NULL)"

Rscript -e "source('$program'); run_mcmc($tid, $nsamp, '$outfile', NULL)"

echo "end:  " `date`