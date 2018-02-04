#!/bin/bash
#$ -V
#$ -N match
#$ -cwd 
#$ -q short-sl7
#$ -l virtual_free=1G,h_rt=0:02:00
#$ -t 1-100
#$ -o ./qsub_out/007-match_libs_$TASK_ID.out.txt
#$ -e ./qsub_out/007-match_libs_$TASK_ID.err.txt

mkdir -p 000-data/007-matched_libs

R CMD BATCH --no-save --no-restore "--args $SGE_TASK_ID" 007-match_libraries.R qsub_out/007-match_libs_$SGE_TASK_ID.Rout.txt

echo "done"