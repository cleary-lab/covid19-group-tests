#! /bin/bash

#$ -cwd
#$ -N Disjunct
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

python disjunct_job.py 18 450 4 1 100 ${SGE_TASK_ID}
