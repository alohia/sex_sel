echo job1; castes=(4 23 5 33 33 44 38 4 14 7 1 44 1);qsub runest.sh ${castes[@]}
echo job2; castes=(23 23 33 4 7 44 44 44 38 38 17 33 2);qsub runest.sh ${castes[@]}
