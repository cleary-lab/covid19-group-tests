#! /bin/bash

#$ -cwd
#$ -N PoolSize
#$ -t 25-36
#$ -q broad
#$ -P regevlab
#$ -l h_vmem=20g
#$ -l h_rt=12:00:00
#$ -e Logs/
#$ -o Logs/

sleep $((SGE_TASK_ID%60))
source /broad/software/scripts/useuse
source /home/unix/bcleary/.my.bashrc
export OMP_NUM_THREADS=1

BASE=../../../
OUT=Results_500k_pop_v2
read dilution m n <<< $(sed -n ${SGE_TASK_ID}p pool_params.txt)
q=$((dilution*m/n))

python -u test_on_simulated_population.py --viral-load-matrix $BASE/Simulated_populations/viral_loads.txt --infection-time $BASE/Simulated_populations/infection_times.txt --onset-time $BASE/Simulated_populations/onset_times.txt --expected-pos-prob 0.95 --savepath $BASE/$OUT/ --n-individuals $n --m-pools $m --max-dilution $dilution | tee $BASE/$OUT/log.n-${n}_m-${m}_q-${q}.txt
