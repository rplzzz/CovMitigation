#!/bin/zsh
#SBATCH -t 2500
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -A clinical_analytics_lab

module load intel/18.0
module load intelmpi/18.0
module load R/3.6.3

echo "start:  " `date`

#program=`Rscript -e 'system.file("scripts/run-mcmc.R", package="CovMitigation")'`
program="./run-mcmc.R"

tid=$SLURM_ARRAY_TASK_ID
nsamp=100000
outfile="mcmc-$tid.rds"
infile="mcmc-$tid.rds"

echo "Run command:"
echo "source('$program'); run_mcmc($tid, $nsamp, '$outfile', '$infile', usescl=TRUE, nproc=$SLURM_CPUS_PER_TASK)"

Rscript -e "source('$program'); run_mcmc($tid, $nsamp, '$outfile', '$infile', usescl=TRUE, nproc=$SLURM_CPUS_PER_TASK)"

echo "end:  " `date`
