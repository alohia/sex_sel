#!/bin/bash

#$ -M akhil.lohia@upf.edu
#$ -m abe
#$ -r n
# -o $HOME/out.txt
# -e $HOME/error.txt
# -pe threaded 4
# -l h_vmem=6G
#$ -cwd

module load matlab

matlab -nodisplay -logfile log.log -r "castes = [$1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11]; disp(castes)"
