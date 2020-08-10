#!/bin/zsh
#SBATCH -t 250
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -A clinical_analytics_lab

module purge
module load intel/18.0
module load intelmpi/18.0
module load R/3.6.3

echo "start:  " `date`

#program=`Rscript -e 'system.file("scripts/run-mcmc.R", package="CovMitigation")'`
program="./mcmc-init.R"

echo "hostname: "
echo `hostname`

tid=$SLURM_ARRAY_TASK_ID
outfile="mcmc-init-$tid.rds"

echo "Run command:"
echo "source('$program'); mcmc_init($tid, '$outfile', $SLURM_CPUS_PER_TASK)"

Rscript -e "source('./mcmc-init.R'); print(system.time(mcmc_init($tid, '$outfile', $SLURM_CPUS_PER_TASK)))"

echo "end:  " `date`
