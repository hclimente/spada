#!/bin/bash
#$ -q normal 
#$ -cwd 
#$ -V 
#$ -wd /data/users/hector/smartas/

source ~/.bashrc
wd="/data/users/hector/smartas"

tumor=`head -n$SGE_TASK_ID $wd/tcga_tags.txt | tail -n1`

python $wd/pipeline/smartas.py -s $JOB_NAME  -t $tumor -wd $wd