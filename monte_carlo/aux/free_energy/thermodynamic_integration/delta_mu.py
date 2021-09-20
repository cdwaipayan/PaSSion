import numpy as np
from scipy.interpolate import BSpline as spline
import matplotlib.pyplot as plt
import argparse

f1 = open('cubic_diamond/N216/P0p03/mu.dat', mode = 'r')                                                                                                          
f2 = open('fluid/N216/P0p03/mu.dat', mode = 'r')   

t = []
mu = []

for line in f1:
    l = line.split()
    t.append( float(l[0]) )
    mu.append( float(l[1]) )

xtal = spline(t, mu, 2)

t = []
mu = []
for line in f2:
    l = line.split()
    t.append( float(l[0]) )
    mu.append( float(l[1]) )

fluid = spline(t, mu, 2)

temps = np.linspace(0.1,0.2,1001)

for ti in temps:
    print(ti,xtal(ti)-fluid(ti))