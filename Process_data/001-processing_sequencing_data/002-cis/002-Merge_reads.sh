#!/bin/bash
#$ -V
#$ -N pearFOS
#$ -cwd
#$ -q short-sl7
#$ -l virtual_free=20G,h_rt=6:00:00
#$ -t 1-6
#$ -o ./qsub_out/002-merging_$TASK_ID.out.txt
#$ -e ./qsub_out/002-merging_$TASK_ID.err.txt

mkdir -p 002-Merged_seq

FORWARD='000-sequences/C8B0WANXX_FOSintra_16s003374-1-1_Li_lane6_'$SGE_TASK_ID'_1_sequence.txt.gz'
REVERSE='000-sequences/C8B0WANXX_FOSintra_16s003374-1-1_Li_lane6_'$SGE_TASK_ID'_2_sequence.txt.gz'
MERGED='002-Merged_seq/FOSintra_lib'$SGE_TASK_ID'_merged.txt'

PEAR=/software/bl/el6.5/PEAR/pear-0.9.6-bin-64/pear-0.9.6-bin-64


$PEAR -f $FORWARD -r $REVERSE -m 150 -v 100 -n 150 -j 1 -p 0.05 -o $MERGED

echo "done"
