#! /bin/bash

#$ -cwd
#$ -N Prevalence
#$ -t 1-32
#$ -q broad
#$ -P regevlab
#$ -l h_vmem=20g
#$ -l h_rt=47:00:00
#$ -e Logs/
#$ -o Logs/

sleep $((SGE_TASK_ID%60))
source /broad/software/scripts/useuse
source /home/unix/bcleary/.my.bashrc
export OMP_NUM_THREADS=1

BASE=../../../
OUT=Results_prevalence
read batch_size num_batches start end <<< $(sed -n ${SGE_TASK_ID}p prevalence_params.txt)

python -u prevalence_model.py --viral-load-matrix $BASE/Simulated_populations/seir_viral_loads.csv --batch-size $batch_size --num-batches $num_batches --start-time $start --end-time $end --num-trials 100 --savepath $BASE/$OUT/ | tee $BASE/$OUT/log.b-${num_batches}_n-${batch_size}_t0-${start}_t1-${end}.txt



