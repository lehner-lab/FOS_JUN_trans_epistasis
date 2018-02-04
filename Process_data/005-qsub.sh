#!/bin/bash
#$ -V
#$ -N CV
#$ -cwd 
#$ -q short-sl7
#$ -l virtual_free=1G,h_rt=0:10:00
#$ -t 1-4000
#$ -o ./qsub_out/005-CV_thermo_$TASK_ID.out.txt
#$ -e ./qsub_out/005-CV_thermo_$TASK_ID.err.txt

mkdir -p 000-data/005-CV_thermo

R CMD BATCH --no-save --no-restore "--args $SGE_TASK_ID" 005-cross_validation_thermo_model.R qsub_out/005-CV_thermo_$SGE_TASK_ID.Rout.txt

echo "done"