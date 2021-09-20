import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import argparse

def main(args):
    f1 = open('rho_temp.dat', mode = 'r')                                                                                                          
    f2 = open('ene_temp.dat', mode = 'r')                                                                                                                                                                                               

    h_nkbt = []
    t      = []
    G      = []

    temp_ref = args.temp
    press    = args.pres
    G_ref    = args.refmu

    count = 0

    for l1, l2 in zip(f1,f2): 
    
        tem = float(l1.split()[0])

        vol = 1./float(l1.split()[1]) 
        ene = float(l2.split()[1]) 
        h   = ( ene + press*vol ) / tem**2 
        h_nkbt.append(h)
        t.append(tem)

    cs = CubicSpline(x=t, y=h_nkbt)
    temps = np.linspace(np.min(t),np.max(t),1001)

    for i,temp in enumerate(temps):
        a  = temp_ref
        b  = temp

        if(b>a):
            G.append( G_ref - cs.integrate(a=a, b=b) )
        elif(b<a):
            G.append( G_ref + cs.integrate(a=b, b=a) )
        else:
            G.append(G_ref)

        print(temp, G[i])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-temp', '--temp', type=float, default=0.0, help='float giving the temperature of the reference system. Default is zero.')
    parser.add_argument('-pres', '--pres', type=float, default=0.0, help='float giving the pressure of the system. Default is zero.')
    parser.add_argument('-refmu', '--refmu', type=float, default=0.0, help='float giving the reference chemical potential. Default is zero.')
    args = parser.parse_args()
    main(args)