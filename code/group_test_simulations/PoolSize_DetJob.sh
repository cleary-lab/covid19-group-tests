#! /bin/bash

#$ -N PoolSizeDet
#$ -l m_mem_free=20G

sleep $((SGE_TASK_ID%60))
conda activate covid-group-tests
export OMP_NUM_THREADS=1

BASE=../../../covid19-group-tests-data

OUT1=Results_det1
DESIGN1="Pooling (higher prevalence, 96 samples, 6 pools).npy"
python -u test_on_simulated_population.py --viral-load-matrix $BASE/Simulated_populations/seir_viral_loads.identification.csv --infection-time $BASE/Simulated_populations/infection_times.csv --onset-time $BASE/Simulated_populations/onset_times.csv --expected-pos-prob 0.95 --savepath $BASE/$OUT1/ --start-time 0 --end-time 120 --pool-compositions "$BASE/$DESIGN1"

OUT2=Results_det2
DESIGN2="Pooling (lower prevalence, 192 samples, 12 pools).npy"
python -u test_on_simulated_population.py --viral-load-matrix $BASE/Simulated_populations/seir_viral_loads.identification.csv --infection-time $BASE/Simulated_populations/infection_times.csv --onset-time $BASE/Simulated_populations/onset_times.csv --expected-pos-prob 0.95 --savepath $BASE/$OUT2/ --start-time 0 --end-time 120 --pool-compositions "$BASE/$DESIGN2"
