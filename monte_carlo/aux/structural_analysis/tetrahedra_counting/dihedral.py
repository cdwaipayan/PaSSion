import numpy as np
from scipy.spatial import cKDTree as KDTree
import boo
import math
from numba import jit
from math import floor
from matplotlib import pyplot as plt
from matplotlib import ticker

# =========================================================================================================================
# =======================  SET UP THE STYLE OF MATPLOTLIB FIGURE TO BE GENERATED ==========================================
# =========================================================================================================================
# Define the style of the plot
# A list of default styles is given at https://matplotlib.org/gallery/style_sheets/style_sheets_reference.html
plt.style.use('bmh')

# Define fonts to be used as standard for the plot
# The Computer Modern roman font (cmr10) is used here to provide a close match to LaTeX text. 
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'cmr10'
plt.rcParams['font.weight'] = 'light'

# Set the width and colour of the border around the plot
plt.rcParams['axes.linewidth'] = 3.5
plt.rc('axes',edgecolor='black')

# Set up figure environment
fig, ax = plt.subplots()

# Set the size of the axes parameters (i.e., the x and y values)
plt.tick_params(labelsize=24)
# Add or remove grid lines from the plot
ax.grid(False)
# Define a vertical and horizontal line
# plt.axvline(x=0.12, linestyle=(0, (1, 1)), color='black', alpha=0.65, linewidth=2.5)
# plt.axhline(y=0.30, linestyle=(0, (1, 1)), color='black', alpha=0.65, linewidth=2.5)

# Set the width of the ticks on the x and y axes
ax.xaxis.set_tick_params(width=3.5, length=5.0, pad=10)
ax.yaxis.set_tick_params(width=3.5, length=5.0)

# Define the title of the plot and the labels for the x and y axes
# plt.title(r'$\bar{q}$ calculations for strong tetrahedra only with $\epsilon_{AA}=1, \epsilon_{BB}=5$', size=26)
ax.set_ylabel(r'$\mathcal{P}(\phi)$', size=34)
ax.set_xlabel(r'$\phi/^{\circ}$', size=34)

# Set the facecolour of the plot to white, as it is grey for the bmh style
ax.set_facecolor('xkcd:white')

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
def periodic_neighbours(pos, ortn, maxdist, L, sphcyl, lstar=0.0):
    maxdistsq = maxdist**2
    rL = 1./L
    bonds = []
    dists = []
    for i in range(len(pos)-1):
        for j in range(i+1, len(pos)):
            distsq = 0
            rij = pos[i,:] - pos[j,:]
            
            for d in range(pos.shape[1]):
                rij[d] -= L * floor(rij[d] * rL + 0.5)
                distsq += rij[d]*rij[d]

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
    return np.array(bonds, np.int64).reshape((len(dists),2)), np.sqrt(np.array(dists, np.float64))

def globalq(qlm,neigh,l,npart):
    """Calculate the global order parameter Q from input qlm and bonds"""
    nqlm = 0
    n = 0
    for i in range(npart):
        nqlm += float(len(neigh[i]))*qlm[i]
        n += float(len(neigh[i]))

    Qlm = nqlm/n
    
    return boo.ql(Qlm)

# def bond_angle(r,nn):


sphcyl = False

maxbondlength = 2.4

#L = 10.6810833#5.29114*2

npart = 64
numConfigurations = 200
calculateEvery = 1

#open The file
posfile = open(file='tet.dat', mode = 'r')
volfile = open(file='../vol.dat', mode = 'r')
# ortnfile = open(file='../ortn.dat', mode = 'r')
angsfile = open(file='angs.dat', mode = 'w')
angfile = open(file='bondang.dat', mode = 'w')
qfile = open(file='q.dat', mode = 'w')

fcounter = 0

pct = 0
pht = 0
pit = 0

gq4 = 0.0
gq6 = 0.0

fcounter = 0

theta_list = []
ave_theta = 0.0

for frame in range(numConfigurations):
    # read in the positions
    npart = int(posfile.readline())

    pos = np.zeros((npart,3))
    ortn = np.zeros((npart,3))

    for part in range(npart):
        pos[part] = [float(x) for x in posfile.readline().split()]
    
    L = float(volfile.readline().split()[1])
    #convert the lines array into something useful
    bonds, dists = periodic_neighbours(pos, ortn, maxbondlength, L, sphcyl, 0.0)
    neighboursofi = getNeighboursOfI(npart=len(pos),bonds=bonds)

    q4m = boo.bonds2qlm(pos, bonds, l=4, periods=L)
    q6m = boo.bonds2qlm(pos, bonds, l=6, periods=L)
    Q4 = globalq(q4m,neighboursofi,l=4,npart=npart)
    Q6 = globalq(q6m,neighboursofi,l=6,npart=npart)

    # print(Q4, Q6)
    count = 0
    ave_theta = 0.0
    qfile.write(str(Q4) + " " + str(Q6) + "\n")

    for j in range(npart):
        rj = pos[j]
        nn = neighboursofi[j]
        for i in range(3):
            if(len(nn) > i):
                ri = pos[nn[i]]
                rij = ri - rj 
                rij = rij - (np.rint(rij / L) * L)
                dij = np.sqrt(np.matmul(rij,rij))
                for k in range(i+1,4):
                    if(len(nn) > k):
                        rk = pos[nn[k]]
                        rjk = rk - rj
                        rjk = rjk - (np.rint(rjk / L) * L)
                        djk = np.sqrt(np.matmul(rjk,rjk))
                        theta_ijk = math.acos( np.matmul(rij,rjk) / (dij*djk) )
                        theta_list.append(theta_ijk*(180/math.pi))
                        #ave_theta += theta_ijk*(180/math.pi)
                        #count += 1
                        angsfile.write(str(theta_ijk*(180/math.pi)) + "\n")
    #angfile.write(str(ave_theta/float(count)) + "\n")

mean_theta = np.mean(theta_list)
std_theta  = np.std(theta_list)
angfile.write(str(mean_theta) + " " + str(std_theta))

ax.xaxis.set_tick_params(width=3.5, length=5.0, pad=10)
ax.yaxis.set_tick_params(width=3.5, length=5.0)
ax.xaxis.set_tick_params(width=2.5, length=3.5, which='minor',  pad=10)
ax.yaxis.set_tick_params(width=2.5, length=3.5, which='minor')

ax.xaxis.set_ticks(np.arange(40.0, 171.0, 10.))
ax.yaxis.set_ticks(np.arange(0.0, 0.51, 0.01))
minor_ticks = np.arange(40.0, 171.0, 5.0)
minor_ticksy = np.arange(0.0, 0.51, 0.005)
ax.set_xticks(minor_ticks, minor=True)
ax.set_yticks(minor_ticksy, minor=True)

plt.hist(theta_list, bins=100, normed=True, alpha=0.5,
         histtype='stepfilled', color='steelblue',
         edgecolor='none');
plt.show()