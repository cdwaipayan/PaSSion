import math
import numpy as np
import scipy.integrate as integrate
import scipy.special as special
import matplotlib.pyplot as plt
import collections

def frequency_count(itt, vals, nr_bins, minn=None, maxx=None):
    ret = []
    if minn == None:
        minn = min(itt)
    if maxx == None:
        maxx = max(itt)
    binsize = (maxx - minn) / float(nr_bins) #man, do I hate int division

    #construct bins
    for x in range(0, nr_bins):
        start = minn + x * binsize
        ret.append([start, start+binsize, 0])

    #assign items to bin
    bins_f = []
    hist_f = []
    for item, val in zip(itt,vals):
        for binn in ret:
            if binn[0] <= item < binn[1]:
                binn[2] += val 

    for binn in ret:  
        bins_f.append(binn[0])
        hist_f.append(binn[2])

    return bins_f, hist_f

f = open('dos_12.dat', mode='r')

t    = []
hist = []

M    = []
rho  = []
ene  = []
prb  = []

for ti in f:
    line = ti.split()
    rho.append( 1000./float(line[0]) )
    ene.append( float(line[1]) )
    prb.append( float(line[2]) )

s = 0.5

for r, e in zip(rho,ene):
    M.append(r + s*e/1000.)

t, hist = frequency_count(M, prb, 50)

t    = np.asarray(t)
hist = np.asarray(hist)

t_bar = np.mean(t)
t_var = np.var(t)

bins = (t-t_bar) / np.sqrt(t_var)
hist = hist / np.max(hist)

plt.plot(bins[0:135875],hist/2.3302548567993915)


# Ising Model
p = lambda m : np.exp( -((m**2/1.1341655**2)-1)**2 * (0.158*(m**2/1.1341655**2)+0.776) )
A, err = integrate.quad(p, a=-3, b=3)
M = np.linspace(-3,3,1000)
p_ising = p(M) / A
plt.plot(M,p_ising)

# plt.show()

f_t = open('t.out', mode='w')
for i, j in zip(hist, bins):
    f_t.write(str(j) + " " + str(i/2.3302548567993915) + '\n')

f_t = open('ising.out', mode='w')
for i, j in zip(p_ising,M):
    f_t.write(str(j) + " " + str(i) + '\n')