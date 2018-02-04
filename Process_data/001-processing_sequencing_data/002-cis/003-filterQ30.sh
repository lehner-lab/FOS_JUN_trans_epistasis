#!/bin/bash
#$ -V
#$ -N filter
#$ -cwd 
#$ -q short-sl7
#$ -l virtual_free=20G,h_rt=2:00:00
#$ -t 1-6
#$ -o ./qsub_out/003-filterQ30.$TASK_ID.out.txt
#$ -e ./qsub_out/003-filterQ30.$TASK_ID.err.txt


mkdir -p 003-filtered

perl 003-filter.pl $SGE_TASK_ID 30

echo "done"
