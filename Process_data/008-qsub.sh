#!/bin/bash
#$ -V
#$ -N fit
#$ -cwd 
#$ -q short-sl7
#$ -l virtual_free=1G,h_rt=2:00:00
#$ -t 1-1000
#$ -o ./qsub_out/008-fit_thermo_on_matched_libs_$TASK_ID.out.txt
#$ -e ./qsub_out/008-fit_thermo_on_matched_libs_$TASK_ID.err.txt

mkdir -p 000-data/008-fit_thermo_on_matched_libraries

R CMD BATCH --no-save --no-restore "--args $SGE_TASK_ID" 008-fit_thermo_on_matched_libraries.R qsub_out/008-fit_thermo_on_matched_libs_$SGE_TASK_ID.Rout.txt

echo "done"