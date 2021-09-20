"""
THIS SCRIPT IS TO BE USED TO CALCULATE THE VALUE OF A DEFINITE INTEGRAL FOR AN UNDEFINED
FUNCTION USING GAUSS-LEGENDRE QUADRATURE.
  - By undefined function we mean something like f(X), where we can evaluate the value of the
    function for a given value of X, but we do not know its functional form.
  - There are three steps to evaluating such an integral:
        1. For a given number of nodes and an integration interval we need to calculate the 
           the values of the nodes (according to some set of rules) at which we are to evaluate  
           the function to be integrated (i.e., we need to find the set {X} for f(X)).
        2. Evaluate f(X) for each value in {X}
        3. Calculate the integral of f(X) using the quadrature rule, and the weights as defined by
           Gauss-Legendre rules.

  - Step 2 is not performed in this script, and so this script is to first to be used to calculate the 
    values {X}, then outside of this script f({X}) is to be evaluated and saved to a file named "function.out"
    which has layout 'X f(X)' listed with X in ASCENDING order. Following this, the script is to be run again 
    in order to evaluate the value of the integral. 
"""
import argparse
import numpy as np
import math
from scipy import special

def trans_ref(npart,lam,temp,rho):

    if(temp==0. or lam==0.):
        return 0.
    return -1.5*(npart-1.0)/npart*np.log(temp*math.pi/lam)-3./(npart*2.)*np.log(npart)+np.log(rho)/npart

def ortn_ref(nsites,lam,temp,headtailt):

    if(temp==0. or lam==0.):
        return 0.
    
    if(nsites==2):
        if(headtailt):
            return -np.log((1.-np.exp(-2.*lam/temp))/(2.*lam/temp))
        else:
            return 0.6921 + np.log(lam/temp)

    elif(nsites==4):
        return 0.453 + 1.5*np.log(lam/temp)

def main(args):

    step   = args.step
    n      = args.n
    lammax = args.lammax

    c = np.exp(3.5)

    a = np.log(c)
    b = np.log(lammax+c)
 
    x, w = special.roots_legendre(n)

    if(step == 1):
        new_x = []

        for i,xi in enumerate(x):
            new_xi = (b-a)*xi/2.0 + (a+b)/2.0
            lami   = np.exp(new_xi) - c
            
            new_x.append(lami)
            print(i+1,lami)

    elif(step == 3):
        integral = 0
        st_err   = 0
        func_tr  = []
        ferr_tr  = []
        func_or  = []
        ferr_or  = []
        rigidt   = False

        intgrnd = open("integrand.dat", "a")

        for i in range(n):
            new_xi = (b-a)*x[i]/2.0 + (a+b)/2.0
            with open("lam"+str(i+1)+"/ave_data.dat", "r") as file:
                func_file = list(file)
                last_line = func_file[-1].split()
                if(last_line[0] == "EIN_TR:" or last_line[0] == "WELL:"):
                    func_tr.append( float(last_line[1]) * np.exp(new_xi) )
                    ferr_tr.append( float(last_line[2]) * np.exp(new_xi) )
                    intgrnd.write(str(new_xi) + " " + str(func_tr[i]) + " " + str(ferr_tr[i]) + "\n")
                else:
                    rigidt = True
                    ferr_or.append( float(last_line[2]) * np.exp(new_xi) )
                    func_or.append( float(last_line[1]) * np.exp(new_xi) )
                    scnd_line = func_file[-2].split()
                    func_tr.append( float(scnd_line[1]) * np.exp(new_xi) )
                    ferr_tr.append( float(scnd_line[2]) * np.exp(new_xi) )
                    intgrnd.write(str(new_xi) + " " + str(func_tr[i]) + " " + str(ferr_tr[i]) + " " + str(func_or[i]) + " " + str(ferr_or[i]) + "\n")

    #   Calculate the value of the integral
        for wi,fi,erri in zip(w,func_tr,ferr_tr):
            integral += wi*fi
            st_err   += (wi*erri)**2

        tr_integral = -(b-a)*integral / 2.0
        tr_err   = (b-a)*math.sqrt(st_err) / 2.0

        err_expo = abs(np.floor(np.log10(np.abs(tr_err))).astype(int))+1

        intgrnd.write("dA2_tr = "+str(round(tr_integral,err_expo))+"("+str(int(10**(err_expo)*round(tr_err,err_expo)))+")\n")

        if(rigidt):
            integral = 0.0
            st_err   = 0.0
            for wi,fi,erri in zip(w,func_or,ferr_or):
                integral += wi*fi
                st_err   += (wi*erri)**2

            or_integral = -(b-a)*integral / 2.0
            or_err   = (b-a)*math.sqrt(st_err) / 2.0

            err_expo = abs(np.floor(np.log10(np.abs(or_err))).astype(int))+1

            intgrnd.write("dA2_or = "+str(round(or_integral,err_expo))+"("+str(int(10**(err_expo)*round(or_err,err_expo)))+")\n")
    
        f0 = open('../step1/ave_data.dat', mode='r')
        f0_file = list(f0)
        last_line = f0_file[-1].split()
        dA1 = float(last_line[1])
        dA1_err = float(last_line[2])

        finp = open('../step1/input.inp', mode='r')
        for line in finp:
            l = line.split()
            if(l[0]=='NPART'):
                npart = float(l[1]) 
            elif(l[0]=='TEMP'):
                temp = float(l[1])
            elif(l[0]=='DENSITY'):
                rho = float(l[1])
            elif(l[0]=='FL_FE'):
                lammax = float(l[2])
                if(rigidt): eta_o  = float(l[3])

        tr_ref = trans_ref(npart,lammax,temp,rho)
        if(rigidt):
            or_ref = ortn_ref(2,lammax*eta_o,temp,False)
            free_energy = tr_ref + or_ref + dA1 + tr_integral + or_integral
            free_error  = np.sqrt(dA1_err**2 + tr_err**2 + or_err**2) 
        else:
            free_energy = tr_ref + dA1 + tr_integral
            free_error  = np.sqrt(dA1_err**2 + tr_err**2) 
        print(free_energy,free_error)
    else:
        raise Exception('step should take a value of 1 or 3')

if __name__ == '__main__':
#   Define parameters which can be read in from the command line
    parser = argparse.ArgumentParser()
    parser.add_argument('-step', '--step', type=int, default=0, help='integer giving the current step (must be 1 or 3)')
    parser.add_argument('-n', '--n', type=int, default=20, help='integer giving the number of nodes for the Gauss-Legendre integration')
    parser.add_argument('-lammax', '--lammax', type=float, default=-1.0, help='lower limit of the integration')
    args = parser.parse_args()
    main(args)
