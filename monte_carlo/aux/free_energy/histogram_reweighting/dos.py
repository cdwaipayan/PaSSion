import numpy as np
import copy

f_hist = open("../../merged_histogram.dat", mode='r')

volume = []
energy = []
counts = []

nc = 3200000. 

for line in f_hist:
    l = line.split()
    volume.append(float(l[0]))
    energy.append(float(l[1]))
    counts.append( np.log(float(l[2])/nc) )

print(np.max(counts))

temps  = [0.105,0.106,0.107,0.108,0.109,0.11,0.112]
press  = [0.005,0.006,0.0065,0.007,0.0075,0.008,0.0085,0.009,0.01]

volume = np.asarray(volume)
energy = np.asarray(energy)#+7500.
counts = np.asarray(counts)
temps  = np.asarray(temps)
press  = np.asarray(press)

v_min = np.min(volume);   v_max = np.max(volume)
e_min = np.min(energy);   e_max = np.max(energy)
b_min = 1./np.max(temps); b_max = 1./np.min(temps)
p_min = np.min(press);    p_max = np.max(press)

ec = (e_max+e_min)/2.
vc = (v_max+v_min)/2.

G     = np.zeros((len(temps),len(press)))
G_Old = np.zeros(len(temps)*len(press))
gg    = 1.

readin = True

if(readin):
    fg = open('../../gibbs.dat',mode='r')
    gibbs = []
    for line in fg:
        l = line.split()
        gibbs.append(float(l[2]))
    for ti, t in enumerate(temps):
        for pi, p in enumerate(press):
            c = ti*len(press)+pi
            G[ti][pi] = gibbs[c]

ln_O = np.zeros(len(energy))

while(gg>1.e-7):
    g_min = np.min(G); g_max = np.max(G)
    G_Old = copy.deepcopy(np.reshape(G,np.shape(G_Old)))
#   Calculate ln(Ω) for each energy and volume pair
    for j,(e,v) in enumerate(zip(energy,volume)):
        oi = 0.
        c1 = -b_max*(e+p_max*v)
        for ti, t in enumerate(temps):
            b = 1./t
            for pi, p in enumerate(press):
                e1 = -b*(e+p*v) - c1
                e2 = b*G[ti][pi]
                oi = oi + np.exp(e1+e2)
        ln_O[j] = counts[j] - (c1+np.log(oi))

#   Calculate the Gibbs free energy using ln(Ω)
    c4 = (np.max(ln_O)+np.min(ln_O))/2.
    for ti, t in enumerate(temps):
        b = 1./t
        for pi, p in enumerate(press):
            gi = 0.
            c3 = -b*(ec + p*vc)
            for j,(o,e,v) in enumerate(zip(ln_O,energy,volume)):
                gi = gi + np.exp( (o-c4) + (-b*(e+p*v)-c3) )
            G[ti][pi] = -t*(c3+c4+np.log(gi))
            print(t, p, G[ti][pi])

    GDiff = np.abs(np.reshape(G,np.shape(G_Old))-G_Old)
    gg    = np.sqrt(np.dot(GDiff,GDiff))
    print(gg)

b    = 1./0.106
pre  = [0.005, 0.006, 0.007, 0.0075, 0.0078, 0.0079, 0.008, 0.00805, 0.0081, 0.0082, 0.0083, 0.0084, 0.0086, 0.009, 0.01]

for p in pre:

    tmp = []
    for j,(o,e,v) in enumerate(zip(ln_O,energy,volume)):
        tmp.append( o-b*(e+p*v) )

#   Calculate probabilities for energy, volume pairs.
#   Density of states shifted by an arbitrary constant not accounted for when 
#   performing self-consistent equations. Shift is required in order to allow
#   for the computation of exponential.
    cnst = -np.min(tmp)
    tmp  = np.exp(np.asarray(tmp)+cnst)

#   Calculate probability distributions for volume
    prt = []
    for i in range(2330,4858):
        count = 0.
        for v, pri in zip(volume,tmp):
            if(v==float(i)):
                count = count + pri
        prt.append(count)

    tmp = []

    if(p==0.01):
        f_o = open('dos_0p0100.dat',mode='w')
    else:
        f_o = open('dos_0p00'+str(int(p*10000))+'.dat',mode='w')

    cnst = np.sum(prt)
    for i,j in enumerate(range(2330,4858)):
        f_o.write( str(1000./j) + " " + str( prt[i]/cnst ) + '\n')

    f_o.close()
