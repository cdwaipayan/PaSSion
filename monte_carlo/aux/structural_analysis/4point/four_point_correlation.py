import numpy as np
import math
import itertools as itool
from math import floor

def getNeighbours(npart, pos, non_cubic, rcut, L):

    rij    = np.zeros( (npart,npart,3) )
    r      = np.zeros( (npart,npart) )
    nneigh = [[] for i in range(npart)]

    for i in range(npart-1):

        ri = pos[i]
        
        for j in range(i+1, npart):
            rj = pos[j]

#       minimum image convention for cubic box
            r_ij = rj - ri
            r_ij = r_ij - (np.rint(r_ij / L) * L)

            rij[i][j] =  r_ij 
            rij[j][i] = -r_ij 

            dist = np.linalg.norm(r_ij)
            
            r[i][j] = dist
            r[j][i] = dist
            
            if dist < rcut:
                nneigh[i].append(j)
                nneigh[j].append(i)
                
    return rij, r, nneigh

# Cutoff used to define bonds
rcut = 1.3

nframes = 1000

nptt = False
non_cubic = False

# npart = 1000
# L     = 15.639424154595982
# posfile = open(file='./testpos.dat', mode = 'r')

# npart = 1000
# L     = 11.696070952851
# posfile = open(file='./diamond/finalpos.dat', mode = 'r')

# npart = 1024
# L     = 11.950412657486
# posfile = open(file='./tetrastack/finalpos.dat', mode = 'r')

# npart   = 2000
# L       = 20.245926274351
# posfile = open(file='./rndm_tetrastack/finaltet.dat', mode = 'r')


# npart   = 898
# rcut    = 2.1
# L       = 20.245926274351
# posfile = open(file='./rndm_tetrastack/finaltet.dat', mode = 'r')

rcut    = 2.6
L       = 17.235477520255
posfile = open(file='./tet.dat', mode = 'r')
tet     = True

if(nptt):
    volfile = open(file='../testbox.dat', mode = 'r')

pi    = math.pi
twopi = 2.*math.pi

for frame in range(nframes):

    if(tet):
        npart = int(posfile.readline())

    pos = np.zeros((npart,3))

#   Read in particles for the current frame
    for part in range(npart):
        pos[part]  = [float(x) for x in posfile.readline().split()]

#   If data is from an NPT simulation read in the box length
    if(nptt):
        L = float(volfile.readline().split()[1])
    
    # Extract the pair separation vectors and nearest neighbours for each particle
    rij, dist, nneigh = getNeighbours(npart, pos, non_cubic, rcut, L)

    sph_coords = []

    for i in range(npart):
        nn_i = nneigh[i]
        nnn  = len(nn_i)
#   skip particles that don't have the required no. of nearest neighbours to form a triplet.
        if nnn >= 2:
#       Generate all possible pairs of numbers between 1 & Number of Nearest Neighbours (NNs).
#       We use this to generate all possible 2 element combinations of neighbours from which a 
#       list of reference coordinate frames can be generated for particle i.
            for nn1 in nn_i:
                for nn2 in nn_i:
                    
                    if(nn1==nn2):
                        continue
            
            #   generate list of particle indices in the cluster of interest
                    c = [nn1,nn2]

            #   ###################################################################
            #       EXTRACT THE AXES DEFINING THE LOCAL COORDINATE SYSTEM WHERE
            #       PARTICLE i IS THE ORIGIN. THIS IS TO BE REPEATED FOR ALL
            #       COMBINATIONS OF NEIGHBOURS.
            #   ###################################################################

                #   vector between particles 1 and 3 that lies in the xz-plane
                    b  = rij[i][c[1]] / dist[i][c[1]]
                
                #   vector between particles 1 and 2 that lies along the z-axis
                    az = rij[i][c[0]] / dist[i][c[0]]
                
                #   vector normal to the plane formed by particles 1, 2 and 3 that 
                #   lies along the y-axis.
                    ay = np.cross(b,az)
                    yhat = np.linalg.norm(ay)

                #   Check for linear chains of particles, as they would be unable to define
                #   a set of orthogonal axes.
                    if(round(yhat,6)==0.):
                        continue
                    
                    ay = ay / yhat

                #   vector that lies along the x-axis
                    ax = np.cross(az,ay)

            #   ###################################################################
            #       USING THE LOCAL COORDINATE FRAME CALCULATE THE DISTRIBUTION OF 
            #       SPHERICAL COORDINATES.
            #   ###################################################################

                    for j in range(npart):
                    #   Cycle over indexes which are the same as the current particle forming
                    #   the origin of the local coordinate system.
                        if(j==i):
                            continue

                    #   Extract the distance between particles i and j
                        r    = dist[i][j]
                    
                        if(r < 4.5):

                            if(j==c[0]):
                                theta = 0.
                                phi   = 0.
                            else:
                            #   Extract the pair separation vector Rij between particles i and j
                                r_ij = rij[i][j]
                                z_ij = np.dot(az,r_ij) 

                            #   Compute the angle between Rij and the z-axis of the local coordinate frame.
                                theta = math.acos( round(z_ij/r, 14) )

                            #   Compute the angle between Rij and the x-axis of the local coordinate frame.
                                if(j==c[1]):
                                    phi = 0.
                                else:
                                #   Projection of Rij onto the xy-plane of the local coordinate frame.
                                    r_xy = r_ij - z_ij*az
                                    rxy  = np.linalg.norm(r_xy)

                                    if(rxy == 0.0):
                                        phi = 0.
                                    else:
                                        rxx  = np.dot(r_xy, ax) / rxy
                                        phi = math.acos( round(rxx, 14) )
                                    #   Check if Rij is pointing along the negative component of the y-axis
                                    #   in the local coordinate frame. If so correct the value of the azimuth
                                    #   for a LH coordinate system.
                                        if(np.dot(r_xy, ay) < 0.):
                                            phi = twopi-phi

                        #   Save the spherical coordinates in a list
                            sph_coords.append([r,theta,phi])

            # break

print('Out of the loop')

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

def near( p, pntList, d0 ):
    cnt=0
    for pj in pntList:
        dist=np.linalg.norm( p - pj )
        if dist < d0:
            cnt += 1 - dist/d0
    return cnt


"""
https://stackoverflow.com/questions/22128909/plotting-the-temperature-distribution-on-a-sphere-with-python
"""

r, theta, phi = np.hsplit(np.asarray(sph_coords), 3)

pointList = []
for theta_i, phi_i, r_i in zip(theta,phi,r):
    ris_t = r_i*np.sin(theta_i)
    x = ris_t*np.cos(phi_i)
    y = ris_t*np.sin(phi_i)
    z = r_i*np.cos(theta_i)
    pointList.append([x,y,z])

pointList = np.reshape(np.asarray(pointList),(len(theta),3))


fig = plt.figure()
ax = fig.add_subplot( 1, 1, 1, projection='3d')

u = np.linspace( 0, 2 * np.pi, 180)
v = np.linspace( 0, np.pi, 90 )

# create the sphere surface
XX = np.max(r)*np.outer( np.cos( u ), np.sin( v ) )
YY = np.max(r)*np.outer( np.sin( u ), np.sin( v ) )
ZZ = np.max(r)*np.outer( np.ones( np.size( u ) ), np.cos( v ) )
print(np.shape(XX))
WW = XX.copy()
for i in range( len( XX ) ):
    for j in range( len( XX[0] ) ):
        x = XX[ i, j ]
        y = YY[ i, j ]
        z = ZZ[ i, j ]
        WW[ i, j ] = near( np.array([x,y,z]), pointList, 0.35 )

WW = WW / np.amax( WW )
myheatmap = WW

ax.plot_surface( XX, YY,  ZZ, cstride=1, rstride=1, facecolors=cm.jet( myheatmap ) )
plt.show() 