#! /bin/bash

#$ -cwd
#$ -N Prevalence
#$ -t 1-896
#$ -q broad
#$ -P regevlab
#$ -l h_vmem=10g
#$ -l h_rt=47:00:00
#$ -e Logs/
#$ -o Logs/

sleep $((SGE_TASK_ID%60))
source /broad/software/scripts/useuse
source /home/unix/bcleary/.my.bashrc
conda activate covid-group-tests
export OMP_NUM_THREADS=1

BASE=../../../
OUT=Results_prevalence_Ct_sputum
read batch_size num_batches start end <<< $(sed -n ${SGE_TASK_ID}p prevalence_params.txt)

PERIOD=growth

python -u prevalence_model.py --viral-load-matrix $BASE/Simulated_populations/sputum/seir_viral_loads_sputum-002.viral_loads.test.npz --kde-path $BASE/Simulated_populations/sputum/KDEs_${PERIOD}/seir_kde --batch-size $batch_size --num-batches $num_batches --start-time $start --end-time $end --num-trials 100 --savepath $BASE/${OUT}_${PERIOD}/ | tee $BASE/${OUT}_${PERIOD}/log.b-${num_batches}_n-${batch_size}_t0-${start}_t1-${end}.txt

PERIOD=decline

python -u prevalence_model.py --viral-load-matrix $BASE/Simulated_populations/sputum/seir_viral_loads_sputum-002.viral_loads.test.npz --kde-path $BASE/Simulated_populations/sputum/KDEs_${PERIOD}/seir_kde --batch-size $batch_size --num-batches $num_batches --start-time $start --end-time $end --num-trials 100 --savepath $BASE/${OUT}_${PERIOD}/ | tee $BASE/${OUT}_${PERIOD}/log.b-${num_batches}_n-${batch_size}_t0-${start}_t1-${end}.txt



