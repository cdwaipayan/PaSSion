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

    refsite = np.zeros((1,3))

    refsite[0]= [ 0., 0.,  -1. ]

    sites = np.zeros((1,3))
    for i in range(1):
        sites[i] = np.matmul(refsite[i],rm)

    return sites


def main(args):
    pos_data = find('pos.dat', '../')

    tet_data = open('../tetrahedra/tet.dat', mode='r')
    tet_idsd = open('../tetrahedra/idtet.dat', mode='r')

    # Loop through data files and plot
    r = np.loadtxt(pos_data[0])
    
    npart  = args.npart 
    ndump  = args.ndump 
    boxl   = args.boxl  
    rcut   = args.rcut  
    nsites = args.ortn
    theta  = args.theta
    kfdel  = np.cos(theta*math.pi/180.)
    nptt   = args.npt


    if(nsites==0):
        ortn = False
    else:
        ortn = True
        ortn_data = find('ortn.dat', '../')
        q = np.loadtxt(ortn_data[0])

    if(nptt):
        box_data = open('../box.dat',mode='r')

    f = open('conflink.dat', 'w')

    for k in range(ndump):

        tet_pos = []
        tet_ids = []

        ntet = int( tet_data.readline() )
        dummy = tet_idsd.readline()
        if(nptt):
            lbox = box_data.readline().split()
            boxl = float(lbox[0])

        for i in range(ntet):
            l1 = tet_data.readline().split()
            l2 = tet_idsd.readline().split()[1:5]
            tet_pos.append([float(l1[0]),float(l1[1]),float(l1[2])])
            pids = []
            for l3 in l2:
                pid = l3.split(',')
                if(len(pid[0].split(']'))==2):
                    pids.append(int(pid[0].split(']')[0]))
                elif(len(pid[0].split('['))==2):
                    pids.append(int(pid[0].split('[')[1]))
                else:
                    pids.append(int(pid[0]))
            tet_ids.append(pids)
        dummy = tet_idsd.readline()
        
        A = np.zeros((ntet,ntet))

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
                    rhat = rij/dist
                    for nsi in range(nsites):
                        earij = -np.dot(si[nsi],rhat)
                        if(earij > kfdel): 
                            for nsj in range(nsites):
                                ebrji = np.dot(sj[nsj],rhat)
                                if(ebrji > kfdel): 
                                    tet_i = 0
                                    tet_j = 0

                                    check1 = False
                                    check2 = False

                                    for ti, tpids in enumerate(tet_ids):
                                        if(i in tpids):
                                            tet_i = ti
                                            check1 = True
                                        elif(j in tpids):
                                            tet_j = ti
                                            check2 = True
                                            break
                                    
                                    if( tet_i != tet_j and (check1 and check2) ):
                                        A[tet_i][tet_j] = 1
                                        A[tet_j][tet_i] = 1

        f.write( "%i, %5.14f\n" % (ntet, boxl) )
        for i in range(ntet):
            f.write( "%5.4f, %5.4f, %5.4f\n" % (tet_pos[i][0],tet_pos[i][1], tet_pos[i][2]) )
        for i in range(ntet):
            neigh = []
            for j in range(ntet):
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
    parser.add_argument('-npart', '--npart', type=int,   default=1,   help='integer giving the number of particles in the system')
    parser.add_argument('-ndump', '--ndump', type=int,   default=1,   help='integer giving the number of simulation snapshots')
    parser.add_argument('-ortn',  '--ortn',  type=int,   default=0,   help='integer giving the rigid body sites for the particle (optional)')
    parser.add_argument('-boxl',  '--boxl',  type=float, default=0.,  help='float giving the cubic box length')
    parser.add_argument('-rcut',  '--rcut',  type=float, default=1.,  help='float giving the radial cutoff for nearest neighbours')
    parser.add_argument('-theta', '--theta', type=float, default=10., help='float giving the patch half-angle for patch B of patchy particles.')
    parser.add_argument('-npt',   '--npt',   action='store_true', default=False, help='Boolean set to true if simulation was in NPT ensemble.')
    args = parser.parse_args()
    main(args)