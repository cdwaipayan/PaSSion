f1 = open('tot_dns.dat', mode='r')
f2 = open('tot_ene.dat', mode='r')
vol = []
ene = []
for line1, line2 in zip(f1,f2):
    vol.append(1000./float(line1))
    ene.append(float(line2))

import math
import numpy as np

ne = 500
nv = 3000

H, volume, energy = np.histogram2d(np.asarray(vol), np.asarray(ene), bins=(nv,ne), range=((2000.,5000.),(-8000.,-7500.)))
f = open('hist.out', mode='w')

for i in range(nv):
    for j in range(ne):
        f.write(str(volume[i]) + " " + str(energy[j]) + " " + str(H[i][j]) + '\n')
