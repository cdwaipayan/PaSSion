from math import pi
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import argparse       
               
def main(args):                                                                                   
    f = open('energy.dat', mode = 'r')                                                                                                                                                                                               

    temp = args.temp
    pres = args.pres
    rho  = args.rho 
    bint = args.binary

    eta  = pi*rho/6.
    f_hs = np.log(rho) - 1. + (4*eta-3*eta**2) / (1.-eta)**2

    if(bint):
        f_hs = f_hs + np.log(0.5)

    b = []
    e = []

    for l in f: 
        b.append( float(l.split()[0]) )
        e.append( float(l.split()[1]) )

    cs = CubicSpline(x=b, y=e)

    for i,beta in enumerate(b):
        print(beta, f_hs + cs.integrate(a=0., b=beta))

    print('CHEMICAL POTENTIAL:', temp, f_hs + cs.integrate(a=0., b=1./temp) + pres/(temp*rho))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-temp', '--temp', type=float, default=0.0, help='float giving the temperature of the system. Default is zero.')
    parser.add_argument('-pres', '--pres', type=float, default=0.0, help='float giving the pressure of the system. Default is zero.')
    parser.add_argument('-rho', '--rho', type=float, default=0.0, help='float giving the density of the system. Default is zero.')
    parser.add_argument('-binary', '--binary', action='store_true', default=False, help="True if considering binary system. Default is False.")
    args = parser.parse_args()
    main(args)