#!/bin/zsh
#SBATCH -t 1000
#SBATCH -n 1
#SBATCH -c 40
#SBATCH -A clinical_analytics_lab

module load gcc/7.1.0
module load R/3.6.1

echo "start:  " `date`

program=`Rscript -e 'system.file("scripts/run-mcmc.R", package="CovMitigation")'`

tid=$SLURM_ARRAY_TASKID
nsamp=100000
outfile="mcmc-warmup-$tid.rds"

echo "Run command:"
echo "source('$program'); run_mcmc($tid, $nsamp, '$outfile', NULL)"

Rscript -e "source('$program'); run_mcmc($tid, $nsamp, '$outfile', NULL)"

echo "end:  " `date`
