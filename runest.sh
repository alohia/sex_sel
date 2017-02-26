#!/bin/bash

#$ -M akhil.lohia@upf.edu
#$ -m abe
#$ -r n
#$ -o $HOME/out.txt
#$ -e $HOME/error.txt
#$ -pe threaded 4
#$ -l h_vmem=5.5G
#$ -l h_rt=02:00:00
#$ -cwd

module load matlab

matlab -nodisplay -logfile log.log -r estimation
