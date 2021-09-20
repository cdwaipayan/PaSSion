from math import pi
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import argparse       
               
def main(args):                                                                                   
    f = open('energy.dat', mode = 'r')                                                                                                                                                                                               

    ref  = args.ref 

    b = []
    e = []

    for l in f: 
        b.append( 1./float(l.split()[0]) )
        e.append( float(l.split()[1]) )

    b = b[::-1]
    e = e[::-1]
    cs = CubicSpline(x=b, y=e)

    for beta in reversed(b):
        print(1./beta, ref + cs.integrate(a=b[-1], b=beta))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-ref', '--ref', type=float, default=0.0, help='float giving the density of the system. Default is zero.')
    args = parser.parse_args()
    main(args)