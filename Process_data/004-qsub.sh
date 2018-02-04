#!/bin/bash
#$ -V
#$ -N param_search
#$ -cwd 
#$ -q short-sl7
#$ -l virtual_free=1G,h_rt=0:10:00
#$ -t 1300-2703
#$ -o ./qsub_out/004-param_search_thermo_$TASK_ID.out.txt
#$ -e ./qsub_out/004-param_search_thermo_$TASK_ID.err.txt

mkdir -p 000-data/004-param_search

R CMD BATCH --no-save --no-restore "--args $SGE_TASK_ID" 004-param_search_thermo_thermo_model.R qsub_out/004-param_search_thermo_$SGE_TASK_ID.Rout.txt

echo "done"