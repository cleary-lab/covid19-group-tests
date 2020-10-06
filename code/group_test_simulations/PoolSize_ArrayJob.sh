#! /bin/bash

#$ -cwd
#$ -N PoolSize
#$ -t 1-90
#$ -q broad
#$ -P regevlab
#$ -l h_vmem=20g
#$ -l h_rt=23:00:00
#$ -e Logs/
#$ -o Logs/

sleep $((SGE_TASK_ID%60))
source /broad/software/scripts/useuse
source /home/unix/bcleary/.my.bashrc
conda activate covid-group-tests
export OMP_NUM_THREADS=1

BASE=../../../
OUT=Results_identification_swab_seir
read m n q <<< $(sed -n ${SGE_TASK_ID}p pool_params.txt)

python -u test_on_simulated_population.py --viral-load-matrix $BASE/Simulated_populations/swab/seir_viral_loads_swab.viral_loads.npz --infection-time $BASE/Simulated_populations/swab/used_pars_swab.csv --expected-pos-prob 0.95 --savepath $BASE/$OUT/ --n-individuals $n --m-pools $m --q-split $q --start-time 10 --end-time 120 | tee $BASE/$OUT/log.n-${n}_m-${m}_q-${q}.txt
