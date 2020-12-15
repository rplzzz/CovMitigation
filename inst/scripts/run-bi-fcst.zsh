#!/bin/zsh
#SBATCH -t 500
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -A clinical_analytics_lab

module load gcc/9.2.0  openmpi/3.1.6
module load python/3.7.7

echo START
date

workdir=`pwd`

cd /home/rl8z/wrk/covid-hosp-sim/python_files
./main_covid_hospitalization_forecast.py --hosp --nc 4 BI-2020-12-08_Adaptive_hospitalized.csv 5000 BI-Adapt-hosp-2020-12-08

./main_covid_hospitalization_forecast.py --hosp --nc 4 BI-2020-12-08_Adaptive-LessControl_hospitalized.csv 5000 BI-AdaptLC-hosp-2020-12-08

./main_covid_hospitalization_forecast.py --hosp --nc 4 BI-2020-12-08_Adaptive-MoreControl_hospitalized.csv 5000 BI-AdaptMC-hosp-2020-12-08

./main_covid_hospitalization_forecast.py --nc 4 BI-2020-12-08_Adaptive_confirmed.csv 5000 BI-Adapt-conf-2020-12-08

./main_covid_hospitalization_forecast.py --nc 4 BI-2020-12-08_Adaptive-LessControl_confirmed.csv 5000 BI-AdaptLC-conf-2020-12-08

./main_covid_hospitalization_forecast.py --nc 4 BI-2020-12-08_Adaptive-MoreControl_confirmed.csv 5000 BI-AdaptMC-conf-2020-12-08


./simulation_aggregation.py $workdir/hospital_census_simulation_results_BI-Adapt-conf-2020-12-08_all.pkl 1 BI-Adapt 0
./simulation_aggregation.py $workdir/hospital_census_simulation_results_BI-AdaptMC-conf-2020-12-08_all.pkl 1 BI-AdaptMC 0
./simulation_aggregation.py $workdir/hospital_census_simulation_results_BI-AdaptLC-conf-2020-12-08_all.pkl 1 BI-AdaptLC 0

echo FIN
date
