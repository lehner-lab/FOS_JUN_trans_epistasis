#!/bin/bash
#$ -V
#$ -N fastqc
#$ -cwd 
#$ -q short-sl7
#$ -l virtual_free=10G,h_rt=1:00:00
#$ -t 1-12
#$ -o ./qsub_out/001-FastQC_$TASK_ID.out.txt
#$ -e ./qsub_out/001-FastQC_$TASK_ID.err.txt

mkdir -p 001-FastQC

INPUT_FILE_LIST=(`ls 000-sequences`)

SGE_TASK_ID_MINUSONE=$((SGE_TASK_ID-1))

/software/bl/el6.5/fastqc/FastQC/fastqc -o 001-FastQC 000-sequences/${INPUT_FILE_LIST[$SGE_TASK_ID_MINUSONE]}

echo "done"
