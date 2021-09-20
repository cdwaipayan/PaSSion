import os
import argparse
import numpy as np
import math
from scipy import special
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
from scipy.interpolate import Akima1DInterpolator as AkimaSpline
from scipy import integrate
from scipy.interpolate import lagrange
import matplotlib.pyplot as plt

nsites    = 2
headtailt = False

def h(w,vv0):

    if(w > 500.):
        val = 3./w - 1.
    else:
        ew = np.exp(w)
        num = 3.*vv0*( 6.+ 2.*ew*(-3.+w) + 4*w + w**2 )
        dem = w*( -w**3 + vv0*( 6. - 6.*ew + 6*w + 3*w**2 + w**3 ) )
        val = num/dem
    
    return val

def ortn_ref(nsites,lam,temp):

    if(temp==0. or lam==0.):
        return 0.
    
    if(nsites==2):
        if(headtailt):
            return -np.log((1.-np.exp(-2.*lam/temp))/(2.*lam/temp))
        else:
            return 0.6921 + np.log(lam/temp)

    elif(nsites==4):
        return 0.453 + 1.5*np.log(lam/temp)

integrand  = []
integrando = []

nlam   = 0
rigidt = False

for i,j,y in os.walk('./step2/'):
    il   = i.split('/')[-1].split('lam0')[-1]
    chil = i.split('/')[-1].split('1')[0]
    if(il=='' or chil=='Lam'):
        continue
    jl = il.split('lam')[-1]
    clam = int(jl)
    if(clam > nlam):
        nlam = clam

lam  = []
Href = []
Oref = []

f0 = open('step1/ave_data.dat', mode='r')
f0_file = list(f0)
last_line = f0_file[-1].split()
dA1 = float(last_line[1])

finp = open('step1/input.inp', mode='r')
for line in finp:
    l = line.split()
    if(l[0]=='NPART'):
        npart = float(l[1]) 
    elif(l[0]=='TEMP'):
        temp = float(l[1])
    elif(l[0]=='DENSITY'):
        rho = float(l[1])
    elif(l[0]=='SCHILLING'):
        rc = float(l[2])
        lammax = float(l[3])
        if(len(l)>4):
            rigidt = True
            if(len(l)>5):
                eta_o  = float(l[-2])
            else:
                eta_o  = float(l[-1])


print("##########################")
print("N  = "+str(int(npart)))
print("T  = "+str(temp))
print("ρ  = "+str(rho))
print("κ  = "+str(lammax))
if(rigidt):
    print("η  = "+str(eta_o))
print("rc = "+str(rc))
print("##########################")

v0    = (4./3.)*math.pi*rc**3
v     = npart / rho
vv0   = v0/v

for i in range(1,nlam+1):
    if(i<10):
        f1 = open('step2/lam00'+str(i)+'/ave_data.dat', mode='r')
        f2 = open('step2/lam00'+str(i)+'/input.inp', mode='r')
    elif(i<100):
        f1 = open('step2/lam0'+str(i)+'/ave_data.dat', mode='r')
        f2 = open('step2/lam0'+str(i)+'/input.inp', mode='r')
    else:
        f1 = open('step2/lam'+str(i)+'/ave_data.dat', mode='r')
        f2 = open('step2/lam'+str(i)+'/input.inp', mode='r')

    f1_file = list(f1)
    last_line = f1_file[-1].split()

    if(rigidt):
        scnd_line = f1_file[-2].split()
        Oref.append( float(last_line[1]) )
        Href.append( float(scnd_line[1]) )
    else:
        Href.append( float(last_line[1]) )

    for line in f2:
        l = line.split()
        if(l[0]=='SCHILLING'):
            lam.append( np.log(float(l[3])) )

count = 0
fi = open('integrand.dat',mode='w')

for lami, Hi in zip(lam,Href):
    href = h(np.exp(lami)/temp, vv0)
    integrand.append( -np.exp(lami)*(Hi-href/temp) )
    fi.write(str(lami) + " " + str(integrand[count]) + "\n")
    count+=1

if(rigidt):
   count = 0
   fio = open('ortn_integrand.dat',mode='w') 
   for lami, Oi in zip(lam,Oref):
        integrando.append( -np.exp(lami)*Oi )
        fio.write(str(lami) + " " + str(integrando[count]) + "\n")
        count+=1
# ####################################################################

# Ideal gas free energy
igas = np.log(rho)-1.
#   Integrate translational component using Akima spline integration
cst  = AkimaSpline(x=lam, y=integrand)
ssex = cst.integrate(a=lam[0], b=lam[-1])
if(rigidt):
#   Compute total free energy
    oref = ortn_ref(nsites,eta_o*lammax,temp)
    cso  = AkimaSpline(x=lam, y=integrando)
    osex = cso.integrate(a=lam[0], b=lam[-1])
    print('SCHILLING-SCHMID FREE ENERGY: ', oref, ssex, osex, dA1, ssex + osex + dA1 + oref + igas)
else:
    print('SCHILLING-SCHMID FREE ENERGY: ', ssex, ssex + dA1, ssex + dA1 + igas  )
    # ####################################################################
    # Calculate free energy using the Carnahan-Starling equation of state
    # ####################################################################
    eta  = math.pi*rho/6.
    hsex = (4.*eta-3*eta**2)/(1.-eta)**2
    print('CARNAHAN-STARLING FREE ENERGY:', hsex, igas+hsex )

# for i in range(1000):
#     r = i/1000.
#     eta  = math.pi*r/6.
#     igas = np.log(r)-1.
#     hsex = (4.*eta-3*eta**2)/(1.-eta)**2
#     print(r, igas+hsex)