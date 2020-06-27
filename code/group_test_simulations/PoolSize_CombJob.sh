#! /bin/bash

#$ -N PoolSizeComb
#$ -t 1-16
#$ -l m_mem_free=20G

sleep $((SGE_TASK_ID%60))
conda activate covid-group-tests
export OMP_NUM_THREADS=1

BASE=../../../covid19-group-tests-data
OUT=Results_comb
read n m q <<< $(sed -n ${SGE_TASK_ID}p comb_params.txt)

python -u test_on_simulated_population.py --viral-load-matrix $BASE/Simulated_populations/seir_viral_loads.identification.csv --infection-time $BASE/Simulated_populations/infection_times.csv --onset-time $BASE/Simulated_populations/onset_times.csv --expected-pos-prob 0.95 --savepath $BASE/$OUT/ --start-time 0 --end-time 120 --pool-compositions $BASE/Design_comb/Acomb.n-${n}_m-${m}_q-${q}.npy
