import numpy as np

f = open('./counts.out', mode='r')

rings = []

r_min = 4
r_max = 9
npart = 1000

for i in range(r_min,r_max+1):
    rings.append([])

for l in f:
    rs = l.split()
    for i,j in enumerate(range(r_min,r_max+1)):
        rings[i].append(float(rs[i]))

for i, r in enumerate(rings):
    print('Number of ' + str(r_min+i) + ': ' + str(np.mean(r)/npart) + ', ' + str(np.std(r)/npart)) 
