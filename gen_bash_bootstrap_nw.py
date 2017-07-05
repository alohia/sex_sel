import numpy as np
castes = [1, 4, 5, 7, 10, 14, 17, 23, 28, 33, 38, 44]
for i in range(1, 101):
    cur_sample = np.append(np.random.choice(castes, 12), i)
    print('echo job' + str(i) + '; castes=(' + str(" ".join(repr(e) for e in cur_sample)) + ');qsub runest.sh ${castes[@]}')
