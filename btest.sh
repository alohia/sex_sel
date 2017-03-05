echo job1; castes=(38 14 28 28 28 28 28 38 1 23 28 1);qsub runest.sh ${castes[@]}
#sleep 1m
echo job2; castes=(28 28 38 28 38 7 28 28 38 44 28 2);qsub runest.sh ${castes[@]}
