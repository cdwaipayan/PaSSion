import numpy as np
import math
import itertools as itool
from math import floor
from math import pi

N       = 4000
rho     = 0.61705
packing = True

if(packing):
    V = N/(rho*6./pi)
else:
    V = N/rho

L = V**(1./3.)

b = 2.*pi / L

rdffile   = open(file='./rdf.dat', mode = 'r')

x  = []
gr = []
for line in rdffile:
    l = line.split()
    x.append(float(l[0]))
    gr.append(float(l[1]))

dr = x[1]-x[0]

for q in x:
    sum = 0.
    for i,r in enumerate(x):
        sum = sum + r *(gr[i]-1.)*np.sin(q*r)
    print(q, 1.+4.*pi*rho*sum*dr/q)