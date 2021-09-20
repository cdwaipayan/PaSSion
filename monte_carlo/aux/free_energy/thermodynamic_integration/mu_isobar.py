import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import argparse

def main(args):
    f1 = open('rho_temp.dat', mode = 'r')                                                                                                          
    f2 = open('ene_temp.dat', mode = 'r')                                                                                                                                                                                               

    G0   = args.refmu
    pres = args.pres

    h_nkbt = []
    t      = []
    G      = []
    G.append(G0)

    for l1, l2 in zip(f1,f2): 
        tem = float(l1.split()[0])
        vol = 1./float(l1.split()[1]) 
        ene = float(l2.split()[1]) 
        h   = ( ene + pres*vol ) / tem**2 
        h_nkbt.append(h)
        t.append(tem)

    cs = CubicSpline(x=t, y=h_nkbt)
    t      = t[::-1]
    h_nkbt = h_nkbt[::-1]

    temps = np.linspace(np.min(t),np.max(t),1001)

    for i,temp in enumerate(temps):
        a  = temp
        b  = t[0]
        G.append( G[0] + cs.integrate(a=a, b=b) )
        print(temp, G[i+1])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-pres', '--pres', type=float, default=0.0, help='float giving the pressure of the system. Default is zero.')
    parser.add_argument('-refmu', '--refmu', type=float, default=0.0, help='float giving the reference chemical potential. Default is zero.')
    args = parser.parse_args()
    main(args)