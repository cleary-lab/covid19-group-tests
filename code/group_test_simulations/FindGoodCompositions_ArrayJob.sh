#! /bin/bash

#$ -cwd
#$ -N FindGoodCompositions
#$ -t 1-1000
#$ -q broad
#$ -P regevlab
#$ -l h_vmem=1g
#$ -l h_rt=4:00:00
#$ -e Logs/
#$ -o Logs/

sleep $((SGE_TASK_ID%60))
source /broad/software/scripts/useuse
source /home/unix/bcleary/.my.bashrc
export OMP_NUM_THREADS=1

python find_good_compositions.py 19 180 0.01 0.9 250 ${SGE_TASK_ID}
