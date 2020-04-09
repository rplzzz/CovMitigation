#!/bin/zsh
#SBATCH -t 2500
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -A clinical_analytics_lab

module load gcc/7.1.0
module load R/3.6.1

echo "start:  " `date`

program=`Rscript -e 'system.file("scripts/run-mcmc.R", package="CovMitigation")'`

tid=$SLURM_ARRAY_TASK_ID
nsamp=100000
outfile="mcmc-warmup-$tid.rds"

echo "Run command:"
echo "source('$program'); run_mcmc($tid, $nsamp, '$outfile', NULL, $SLURM_CPUS_PER_TASK)"

Rscript -e "source('$program'); run_mcmc($tid, $nsamp, '$outfile', NULL, $SLURM_CPUS_PER_TASK)"

echo "end:  " `date`
