import math
import sys
import os, fnmatch
import numpy as np
import itertools

import argparse

# =========================================================================================================================
# =================================  MAIN BODY OF SCRIPT USED TO PLOT THE DATA ============================================
# =========================================================================================================================

# Method to find files with a given regex defined by pattern in a given directory defined by path
# Can be used to read in multiple without having to explicitly define the filenames
def find(pattern, path):
    result = []
    for root, _, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name)) 
    return result

def q_to_sites ( nsites, q ):
    norm = sum ( q**2 ) # quaternion squared length
    if ( abs ( norm - 1. ) > 1.e-06 ):
        exit("ERROR IN QUATERNION NORM")

    rm = np.zeros((3,3))
    # write out column by column, for clarity
    rm[0] = [ q[0]**2+q[1]**2-q[2]**2-q[3]**2,   2.*( q[1]*q[2]+q[0]*q[3] ),     2.*( q[1]*q[3]-q[0]*q[2] )     ]
    rm[1] = [  2.*( q[1]*q[2]-q[0]*q[3] ),    q[0]**2-q[1]**2+q[2]**2-q[3]**2,   2.*( q[2]*q[3]+q[0]*q[1] )     ]
    rm[2] = [  2.*( q[1]*q[3]+q[0]*q[2] ),     2.*( q[2]*q[3]-q[0]*q[1] ),    q[0]**2-q[1]**2-q[2]**2+q[3]**2   ]

    refsite = np.zeros((nsites,3))

    if(nsites==1):
        refsite[0]= [ 0., 0.,  1. ]

    elif(nsites==2):
        refsite[0]= [ 0., 0., -1. ]
        refsite[1]= [ 0., 0.,  1. ]

    elif(nsites==3):
        refsite[0]= [ 0.0,          0.0,          1.0 ]
        refsite[1]= [ 0.0,   (np.sqrt(3.0)/2.0), -0.5 ]
        refsite[2]= [ 0.0,  -(np.sqrt(3.0)/2.0), -0.5 ]

    elif(nsites==4):
        refsite[0]= [  1./np.sqrt(3.),  1./np.sqrt(3.),  1./np.sqrt(3.) ]
        refsite[1]= [ -1./np.sqrt(3.), -1./np.sqrt(3.),  1./np.sqrt(3.) ]
        refsite[2]= [  1./np.sqrt(3.), -1./np.sqrt(3.), -1./np.sqrt(3.) ]
        refsite[3]= [ -1./np.sqrt(3.),  1./np.sqrt(3.), -1./np.sqrt(3.) ]

    sites = np.zeros((nsites,3))
    for i in range(nsites):
        sites[i] = np.matmul(refsite[i],rm)

    return sites


def main(args):
    pos_data = find('pos.dat', '../')

    # Loop through data files and plot
    r = np.loadtxt(pos_data[0])
    
    npart  = args.npart 
    ndump  = args.ndump 
    boxl   = args.boxl  
    rcut   = args.rcut  
    nsites = args.ortn
    thetaA = args.thetaA
    thetaB = args.thetaB
    kfdelA = np.cos(thetaA*math.pi/180.)
    if(nsites>1 and thetaB > 0.):
        kfdelB = np.cos(thetaB*math.pi/180.)
        kfdel  = [kfdelA,kfdelB]
    else:
        kfdel = []
        for i in range(nsites):
            kfdel.append(kfdelA)


    if(nsites==0):
        ortn = False
    else:
        ortn = True
        ortn_data = find('ortn.dat', '../')
        q = np.loadtxt(ortn_data[0])

    f = open('conflink.dat', 'w')

    for k in range(ndump):
        A = np.zeros((npart,npart))
        for i in range(npart-1):
            j1 = npart*k + i
            ri = r[j1]
            if(ortn):
                si = q_to_sites(nsites,q[j1])
            for j in range(i+1, npart):
                j2 = npart*k + j
                rj = r[j2]
                if(ortn):
                    sj = q_to_sites(nsites,q[j2])
                rij = ri - rj
                rij = rij - (np.rint(rij / boxl) * boxl)
                dist = math.sqrt(np.matmul(rij, rij))

                if(dist <= rcut):
                    if(ortn):
                        rhat = rij/dist
                        for nsi in range(nsites):
                            earij = -np.dot(si[nsi],rhat)
                            if(earij > kfdel[nsi]): 
                                for nsj in range(nsites):
                                    ebrji = np.dot(sj[nsj],rhat)
                                    if(ebrji > kfdel[nsj]): 
                                        A[i][j] = 1
                                        A[j][i] = 1
                    else:
                        A[i][j] = 1
                        A[j][i] = 1

        f.write( "%i, %5.14f\n" % (npart, boxl) )
        for i in range(npart):
            j1 = npart*k + i
            f.write( "%5.4f, %5.4f, %5.4f\n" % (r[j1][0],r[j1][1], r[j1][2]) )
        for i in range(npart):
            neigh = []
            for j in range(npart):
                if(A[i,j]==1):
                    neigh.append(j+1)

            f.write("%i, %i\n" % (i+1, len(neigh)))
            for ney in neigh:
                f.write("%i, " % (ney))
            f.write("\n")

# =============================================================
# Define parameters which can be read in from the command line
# =============================================================
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
# Paramters defining the diamond system
    parser.add_argument('-npart',  '--npart',  type=int,   default=1,  help='integer giving the number of particles in the system')
    parser.add_argument('-ndump',  '--ndump',  type=int,   default=1,  help='integer giving the number of simulation snapshots')
    parser.add_argument('-ortn',   '--ortn',   type=int,   default=0,  help='integer giving the rigid body sites for the particle (optional)')
    parser.add_argument('-boxl',   '--boxl',   type=float, default=0., help='float giving the cubic box length')
    parser.add_argument('-rcut',   '--rcut',   type=float, default=1., help='float giving the radial cutoff for nearest neighbours')
    parser.add_argument('-thetaA', '--thetaA', type=float, default=10.,help='float giving the patch half-angle for patch A of patchy particles')
    parser.add_argument('-thetaB', '--thetaB', type=float, default=-1., help='float giving the patch half-angle for patch B of patchy particles. If not given assumed two equal patches used.')
    args = parser.parse_args()
    main(args)