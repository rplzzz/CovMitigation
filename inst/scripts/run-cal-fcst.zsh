#!/bin/zsh
#SBATCH -t 100
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -A clinical_analytics_lab

module load gcc/9.2.0  openmpi/3.1.6
module load python/3.7.7

echo START
date

workdir=`pwd`
jobtag=2020-12-06

cd /home/rl8z/wrk/covid-hosp-sim/python_files
./main_covid_hospitalization_forecast.py --nc 4 cal-infect-stream-2020-12-06.csv 5000 CAL-2020-12-06


./simulation_aggregation.py $workdir/hospital_census_simulation_results_CAL-2020-12-06_BAU.pkl 1 CAL-BAU-$jobtag 0
./simulation_aggregation.py $workdir/hospital_census_simulation_results_CAL-2020-12-06_Flat_20pct_Surge.pkl 1 CAL-Flat20-$jobtag 0
./simulation_aggregation.py $workdir/hospital_census_simulation_results_CAL-2020-12-06_Seasonality_adjustment.pkl 1 CAL-Season-$jobtag 0
./simulation_aggregation.py $workdir/hospital_census_simulation_results_CAL-2020-12-06_Seasonality_and_Travel.pkl 1 CAL-SeasonTravel-$jobtag 0
./simulation_aggregation.py $workdir/hospital_census_simulation_results_CAL-2020-12-06_Travel_adjustment.pkl 1 CAL-Travel-$jobtag 0




echo FIN
date
