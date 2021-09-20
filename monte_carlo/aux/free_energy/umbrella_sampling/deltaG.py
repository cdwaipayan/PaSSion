import copy
import numpy as np
import matplotlib.pyplot as plt

def diff(n1,g1,n2,g2):
    
    d  = 0.
    gd = 0.
    b  = 0.

    check = True

    counter = 0

    while(check):
        di = 0
        d = 0
        gd = 0
        for ni,gi in zip(n1,g1):
            for nj,gj in zip(n2,g2):
                if(ni == nj):
                    di = gi - gj + b
                    d  += di**2
                    gd += 2.*di

        bo = b
        b  -= 0.01*gd

        gdiff = np.sqrt( (b-bo)**2 )
        counter += 1
        
        if(gdiff<1.e-14 or counter > 100000):
            # print(gd)
            check = False

    return b

finaldg = []
finalnc = []
finaler = []

for i in range(0,261,10):

    # if((i-5)%2==0 and i!=65):
    #     continue

    if(i==0):
        f = open('nmax000/umbrella.dat',mode='r')
    elif(i==5):
        f = open('nmax005/umbrella.dat',mode='r')
    elif(i<100):
        f = open('nmax0'+str(i)+'/umbrella.dat',mode='r')
    else:
        f = open('nmax'+str(i)+'/umbrella.dat',mode='r')

    nc = []
    dg = []
    er = []

    # read = False
    
    for line in reversed(f.readlines()):
        
        if(line.split()[0] == 'STEP:'):
            break
        
        nci = float(line.split()[0])
        p   = float(line.split()[2])
        eri = float(line.split()[3])

        if(p>0.01 and nci > 3.):
            if(i==0):
                dgi = -np.log(p/3375.)
            else:
                dgi = -0.075*(nci-i)**2 - np.log(p/3375.)
            
            nc.append(nci)
            dg.append(dgi)
            er.append(eri/p)

            # print(nci,dgi,eri/p)


    if(i>0):
        b = diff(oldnc,olddg,nc,dg)
        # print(i,oldnc,olddg,nc,dg,b)
        oldnc = nc
        olddg = dg - b
    else:
        oldnc = nc
        olddg = dg

    for f1, f2, f3 in (zip(list(reversed(oldnc)),list(reversed(olddg)),list(reversed(er)))):
        finalnc.append(f1)
        finaldg.append(f2)
        finaler.append(f3)

    plt.errorbar(nc,dg,yerr=er, marker='o',linewidth=1.5,markeredgecolor='black', capsize=3.0)

plt.show()

plt.errorbar(finalnc,finaldg,yerr=finaler, marker='o',linewidth=1.5,markeredgecolor='black', capsize=3.0)

plt.show()

# SORT DATA FIRST
zipped_lists = zip(finalnc,finaldg,finaler)
sorted_lists = sorted(zipped_lists)
tuples = zip(*sorted_lists)
finalnc,finaldg,finaler = [list(tup) for tup in tuples]

for f1, f2, f3 in (zip(finalnc,finaldg,finaler)):
    print(f1,f2,f3)