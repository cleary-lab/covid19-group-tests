#! /bin/bash

#$ -N Design_PBEST
#$ -l m_mem_free=20G

conda activate covid-group-tests
export OMP_NUM_THREADS=1

BASE=../../../../2020.03.group-testing/covid19-group-tests-data
OUT=Results_PBEST

python -u test_on_simulated_population.py --viral-load-matrix $BASE/Simulated_populations/seir_viral_loads.identification.csv --infection-time $BASE/Simulated_populations/infection_times.csv --onset-time $BASE/Simulated_populations/onset_times.csv --expected-pos-prob 0.95 --savepath $BASE/$OUT1/ --start-time 0 --end-time 120 --pool-compositions $BASE/Design_PBEST/poolingMatrix.npy
