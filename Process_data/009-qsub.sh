#!/bin/bash
#$ -V
#$ -N fdr
#$ -cwd 
#$ -q short-sl7
#$ -l virtual_free=1G,h_rt=0:40:00
#$ -t 501-1000
#$ -o ./qsub_out/009-FDR_matched_libs_$TASK_ID.out.txt
#$ -e ./qsub_out/009-FDR_matched_libs_$TASK_ID.err.txt

mkdir -p 000-data/009-fdr_matched_libs

R CMD BATCH --no-save --no-restore "--args $SGE_TASK_ID" 009-compute_thermo_epi_matched_libs.R qsub_out/009-FDR_matched_libs_$SGE_TASK_ID.Rout.txt

echo "done"