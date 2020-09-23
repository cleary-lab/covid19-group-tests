#! /bin/bash

#$ -N Design_factors
#$ -l m_mem_free=20G

conda activate covid-group-tests
export OMP_NUM_THREADS=1

BASE=../../../../2020.03.group-testing/covid19-group-tests-data

OUT=Results_ordpair
python -u test_on_simulated_population.py --viral-load-matrix $BASE/Simulated_populations/seir_viral_loads.identification.csv --infection-time $BASE/Simulated_populations/infection_times.csv --onset-time $BASE/Simulated_populations/onset_times.csv --expected-pos-prob 0.95 --savepath $BASE/$OUT/ --ordered-pairs --n-individuals 96 --m-pools 20 --q-split 2 --start-time 0 --end-time 120
python -u test_on_simulated_population.py --viral-load-matrix $BASE/Simulated_populations/seir_viral_loads.identification.csv --infection-time $BASE/Simulated_populations/infection_times.csv --onset-time $BASE/Simulated_populations/onset_times.csv --expected-pos-prob 0.95 --savepath $BASE/$OUT/ --ordered-pairs --n-individuals 384 --m-pools 40 --q-split 2 --start-time 0 --end-time 120

OUT=Results_3factor
python -u test_on_simulated_population.py --viral-load-matrix $BASE/Simulated_populations/seir_viral_loads.identification.csv --infection-time $BASE/Simulated_populations/infection_times.csv --onset-time $BASE/Simulated_populations/onset_times.csv --expected-pos-prob 0.95 --savepath $BASE/$OUT/ --start-time 0 --end-time 120 --pool-compositions $BASE/Design_3factor/A3factor.n-384_m-42_q-3.npy
python -u test_on_simulated_population.py --viral-load-matrix $BASE/Simulated_populations/seir_viral_loads.identification.csv --infection-time $BASE/Simulated_populations/infection_times.csv --onset-time $BASE/Simulated_populations/onset_times.csv --expected-pos-prob 0.95 --savepath $BASE/$OUT/ --start-time 0 --end-time 120 --pool-compositions $BASE/Design_3factor/A3factor.n-288_m-30_q-3.npy
python -u test_on_simulated_population.py --viral-load-matrix $BASE/Simulated_populations/seir_viral_loads.identification.csv --infection-time $BASE/Simulated_populations/infection_times.csv --onset-time $BASE/Simulated_populations/onset_times.csv --expected-pos-prob 0.95 --savepath $BASE/$OUT/ --start-time 0 --end-time 120 --pool-compositions $BASE/Design_3factor/A3factor.n-96_m-18_q-3.npy
python -u test_on_simulated_population.py --viral-load-matrix $BASE/Simulated_populations/seir_viral_loads.identification.csv --infection-time $BASE/Simulated_populations/infection_times.csv --onset-time $BASE/Simulated_populations/onset_times.csv --expected-pos-prob 0.95 --savepath $BASE/$OUT/ --start-time 0 --end-time 120 --pool-compositions $BASE/Design_3factor/A3factor.n-144_m-30_q-3.npy
