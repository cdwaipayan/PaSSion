import numpy as np
from scipy.spatial import cKDTree as KDTree
import boo
import math
import itertools as itool
from numba import jit
from math import floor
from numpy.linalg import inv

def init_list_of_objects(size):
    list_of_objects = list()
    for i in range(0,size):
        list_of_objects.append( list() ) #different object reference each time
    return list_of_objects

# ! =================================================================================================================
# ! Convert angles in degrees to radians
# ! =================================================================================================================
def to_radians(angles):
    angrad = []
    for i in angles:
        angrad.append(i*np.pi/180.)
    return angrad

def shortd(rij, ui, uj, l_star):

    lami = 0.0
    lamj = 0.0
    
    uri  = np.dot(ui,rij)
    urj  = np.dot(uj,rij)
    uij  = np.dot(ui,uj)
    l2   = l_star/2.0

    cc = (1.0-uij**2)

    if(cc < 1e-06):
        if(uri != 0.0):
            
            lami = np.sign(uri)*l2
            lamj = lami*uij - urj
            
            if(abs(lamj) > l2):
                lamj = np.sign(lamj)*l2

        else:
            lami = 0.0
            lamj = 0.0
        
    else:
        lami = (uri - uij*urj) / cc
        lamj = (uij*uri - urj) / cc

        li = abs(lami) - l2
        lj = abs(lamj) - l2
        
        if(li > 0.0 or lj > 0.0):
            if(li > lj):
                lami = np.sign(lami)*l2
                lamj = lami*uij - urj

                if(abs(lamj) > l2):
                    lamj = np.sign(lamj)*l2

            else:

                lamj = np.sign(lamj)*l2
                lami = lamj*uij + uri

                if(abs(lami) > l2):
                    lami = np.sign(lami)*l2

    
    dij = rij + lamj*uj - lami*ui
    return dij

def quat_to_unit(q):
    u = np.zeros(3)
    u[0] = 2.0*( q[1]*q[3] + q[0]*q[2] )
    u[1] = 2.0*( q[2]*q[3] - q[0]*q[1] )
    u[2] = q[3]**2 + q[0]**2 - q[1]**2 - q[2]**2
    return u

#  =================================================================================================================
#  Convert angle-axis vector into a unit vector using formula for the rotation matrix
#  =================================================================================================================
def aatouvec(v):
    v = np.asarray(v)
    theta = np.sqrt(np.dot(v,v))
    q = np.zeros(4)
    u = np.zeros(3)

    if(theta <= 1.e-12):
        q = [1, 0, 0, 0]
    else:
        thetah = 0.5*theta
        st = np.sin(thetah)
        q[0] = np.cos(thetah)
        q[1:4] = v[:]*st / theta

    u = quat_to_unit(q)

    return u

#  =================================================================================================================
#  Given a set of box parameters: boxp = (a, b, c, α, β, γ) calculate the cell matrix L, and inverse cell matrix U
#  as defined in Equations (2) and (6) of K. D. Gibson, H. A. Scheraga, J Phys Chem, 1995, 99, 3752-3764, 
#  Crystal Packing without Symmetry Constraints. 1.
#  =================================================================================================================
def unit_cell(boxp):

    L = np.zeros((3,3))
    U = np.zeros((3,3))

#  Extract cell lengths from box parameters
    a   = boxp[0] 
    b   = boxp[1] 
    c   = boxp[2]

#  Extract trig values using cell angles from box parameters    
    cosalp = np.cos(boxp[3]); cosbet = np.cos(boxp[4])
    cosgam = np.cos(boxp[5]); singam = np.sin(boxp[5])

#  Unit volume
    W = np.sqrt( 1.0 - cosalp**2 - cosbet**2 - cosgam**2 + (2.0*cosalp*cosbet*cosgam) )

    c3 = c*( (cosalp-cosbet*cosgam) / singam )

#  Calculate upper triangular cell matrix, L
    L[0]  = (a,   b*cosgam, c*cosbet    )
    L[1]  = (0.0, b*singam, c3          )
    L[2]  = (0.0, 0.0     , c*W / singam)

#  Using Eqn. (8) from Crystal Packing without Symmetry Constraints. 1., calculate the inverse cell matrix U
    U = np.transpose( inv(L) )

    return L, U

"""Give a sorted list of neighbours from the bonds calculated by pyboo."""
def getNeighboursOfI(npart,bonds):
    neighboursofi = []
    for part in range(npart):
        ilist = []
        for bond in bonds:
            if(bond[0] == part):
                ilist.append(bond[1])
            elif (bond[1] == part):
                ilist.append(bond[0])
        ilist.sort()
        neighboursofi.append(ilist)
    return neighboursofi

"""Used to calculate neighbours from a cubic box. Copied from the pyboo -glass in periodic boundary conditions- case."""
# @jit(nopython=True) #numba to speed up our brute force approach
def periodic_neighbours(pos, ortn, non_cubic, maxdist, L, sphcyl, lstar=0.0):
    maxdistsq = maxdist**2
    bonds = []
    dists = []
    vecs  = init_list_of_objects(len(pos))

    for i in range(len(pos)-1):
        
        ri = pos[i]
        if(non_cubic): 
            si = np.matmul(inv(L),ri)
        
        for j in range(i+1, len(pos)):
            distsq = 0
            rj = pos[j]
#       minimum image convention for arbitrary unit cell      
            if(non_cubic):
                sj  = np.matmul(inv(L),rj)
                sij = si - sj
                sij -= np.rint(sij)
                rij = np.matmul(L,sij)
            else:
#       minimum image convention for cubic box
                rij = ri - rj
                rij = rij - (np.rint(rij / L) * L)

            distsq = np.matmul(rij,rij)

            if(sphcyl):
                ui = ortn[i,:]
                uj = ortn[j,:]
                diff = shortd(rij, ui, uj, lstar)
                diffsq = np.dot(diff, diff)
            else:
                diffsq = distsq
            
            if diffsq < maxdistsq:
                bonds.append(i)
                bonds.append(j)
                dists.append(distsq)
                vecs[i].append(rij)
                vecs[j].append(-rij)
    return np.array(bonds, np.int64).reshape((len(dists),2)), np.sqrt(np.array(dists, np.float64)), vecs

def globalq(qlm,neigh,l,npart):
    """Calculate the global order parameter Q from input qlm and bonds"""
    nqlm = 0
    n = 0
    for i in range(npart):
        nqlm += float(len(neigh[i]))*qlm[i]
        n += float(len(neigh[i]))

    Qlm = nqlm/n
    
    return boo.ql(Qlm)

# Set sphcyl to true to account for cylindrical symmetry of patchy particles
# lstar defines the length of the effective line segment of the particle
sphcyl = False
lstar  = 0.5

# Cutoff used to define bonds
rcut = 1.3

cm = 16.069202332988

# Threshold with which to define a cluster as being a tetrahedron
qmin = 0.99

nframes = 20000#6000

# Set to true if the rotational coordinates are to be given as unit vectors
unit_vector = False

npart = 1000

xyz = False

posfile   = open(file='../pos.dat', mode = 'r')
ortnfile  = open(file='../ortn.dat', mode = 'r')

lrgst = True
lrgstfile = open(file='../lrgst_x_clstr.xyz',mode='r')

non_cubic = False

viewtet = open(file='./tet.xyz', mode = 'w+')
tet = open(file='./tet.dat', mode = 'w+')
ntet = open(file='./numtet.dat', mode = 'w+')
tetid = open(file='./idtet.dat', mode = 'w+')

tetsum = 0

for frame in range(nframes):

    if(lrgst):
        pos2 = np.zeros((2*npart,3))
        ortn2 = np.zeros((2*npart,3))
        count = 0
        for i in range(3*npart+2):
            part = lrgstfile.readline().split()
            if(i>1):
                if(part[0]=='Li' or part[0]=='He' or part[0]=='H'):
                    pos2[count] = [float(part[1]),float(part[2]),float(part[3])]
                    count = count + 1
    else:
        pos = np.zeros((npart,3))
        ortn = np.zeros((npart,3))

        for part in range(npart):
            pos[part]  = [float(x) for x in posfile.readline().split()]
            q = [float(x) for x in ortnfile.readline().split()]
            if(unit_vector):
                ortn[part] = q[part]
            else:
                if(len(q)==4):
                    ortn[part] = quat_to_unit(q)
                else:
                    ortn[part] = aatouvec(q)
        
        pos2  = []
        ortn2 = []

        for i in range(npart):
            pos2.append( pos[i] + ortn[i]*(lstar/2.0) )
            pos2.append( pos[i] - ortn[i]*(lstar/2.0) )
            ortn2.append(ortn[i])
            ortn2.append(-ortn[i])

        pos2 = np.asarray(pos2)
        ortn2 = np.asarray(ortn2)

    #convert the lines array into something useful
    bonds, dists, rij = periodic_neighbours(pos2, ortn2, non_cubic, rcut, cm, sphcyl, lstar)
    nneigh = getNeighboursOfI(npart=npart*2, bonds=bonds)

    tetrahedra = []

    for i in range(npart*2):
#   skip particles that don't have the required no. of neighbours to form a tetrahedron.
        if len(nneigh[i]) >= 3:
#   generate a list of the positions of nearest neighbours. 
#   This must be done in order to account for periodic boundaries later.
            min_neigh = []
            ri = pos2[i]
            min_neigh.append(ri)
            
#       Extract position of nearest neighbours using the minimum image convention
            for j in range(len(nneigh[i])):
                min_neigh.append(ri - rij[i][j])

#       generate all possible 3 element combinations of the neighbours
            for j in itool.combinations(range(1,len(min_neigh)), 3):
#           generate list of particle indices in the cluster of interest
                c   = [0, j[0], j[1], j[2]]
                clu = [ i, nneigh[i][j[0]-1], nneigh[i][j[1]-1], nneigh[i][j[2]-1] ]
#           calculate the center of mass of the cluster.
#           this corresponds to a candidate tetrahedral center.
                com = 0.25*(min_neigh[c[0]] + min_neigh[c[1]] + min_neigh[c[2]] + min_neigh[c[3]] )

#           calculate tetrahedral order parameter q
                q = 0
                for k in range(3):
                    rk = min_neigh[c[k]]
                    rkc = com - rk
                    dkc = np.sqrt(np.matmul(rkc,rkc))

                    for l in range(k+1,4):
                        rl = min_neigh[c[l]]
                        rcl = com - rl
                        dcl = np.sqrt(np.matmul(rcl,rcl))
                        
                        cos_phi_kl = np.matmul(rkc,rcl) / (dkc*dcl)
                        q += ( cos_phi_kl + (1.0/3.0) )**2

                q = 1.0 - (3.0/8.0)*q

#           if the calculated q satisfies the minimum threshold, save the current cluster as a tetrahedron
                if( q >= qmin ):
                    tetrahedra.append([q,clu,com])

#           Remove particle i from the nearest neighbour list of all of its neighbours.
#           This doubly serves the purpose of saving computional time, and prevents overcounting tetrahedra.
            for j in nneigh[i]:
                index = nneigh[j].index(i)
                nneigh[j].remove(i)
                del(rij[j][index])

    tetid.write("Configuration "+str(frame)+ ", No. Tetrahedra: " + str(len(tetrahedra)) +"\n")
    for itetr in tetrahedra:
        q = itetr[0]
        c = itetr[1]
        tetid.write(str(q) + " " + str(c) + "\n")
    tetid.write("\n")
    
    tet.write(str(len(tetrahedra))+"\n")
    for itetr in tetrahedra:
        ri = itetr[2]
        tet.write(str(ri[0]) + " " + str(ri[1]) + " " + str(ri[2])+"\n")
    
    ntet.write(str(len(tetrahedra))+"\n")
    tetsum += len(tetrahedra)

    if(xyz and frame%500==0):
        viewtet.write(str(len(pos)*3+len(tetrahedra))+"\n")
        viewtet.write("\n")
        for ri, ui in zip(pos,ortn):
            viewtet.write( "Au " + str(ri[0]) + " " + str(ri[1]) + " " + str(ri[2]) + "\n")
            viewtet.write( "O " + str(ri[0] + 0.5*lstar*ui[0]) + " " + str(ri[1] + 0.5*lstar*ui[1]) + " " + str(ri[2] + 0.5*lstar*ui[2]) + "\n")
            viewtet.write( "N " + str(ri[0] - 0.5*lstar*ui[0]) + " " + str(ri[1] - 0.5*lstar*ui[1]) + " " + str(ri[2] - 0.5*lstar*ui[2]) + "\n")

        for i in tetrahedra:
            ri = i[2]
            viewtet.write( "C " + str(ri[0]) + " " + str(ri[1]) + " " + str(ri[2]) + "\n")

ntet.write("Average Number of Tetrahedra: " + str(float(tetsum)/float(nframes))+"\n")
