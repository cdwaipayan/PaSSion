import numpy as np
from scipy.spatial import cKDTree as KDTree
import boo
import math
import itertools as itool
from numba import jit
from math import floor
from numpy.linalg import inv

import matplotlib.pyplot as plt
from matplotlib import font_manager as fm, rcParams
import matplotlib.ticker as ticker

import argparse

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
ax.set_xlabel(r'$r / \sigma$', size=34)
ax.set_ylabel(r'$g(r)$', size=34, rotation=90)

# Set the facecolour of the plot to white, as it is grey for the bmh style
ax.set_facecolor('xkcd:white')
# Define a list of colours to be used for the points on the scatter plot
# List of colours taken from https://matplotlib.org/users/colors.html (can be expanded if required)
colour = 0
colours = ['teal', 'coral', 'lavender', 'aqua', 'olive', 'gold']
# =========================================================================================================================
# =========================================================================================================================

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
    L[0]  = (a,   0.0, 0.0    )
    L[1]  = (b*cosgam, b*singam, 0.0          )
    L[2]  = (c*cosbet, c3     , c*W / singam)

#  Using Eqn. (8) from Crystal Packing without Symmetry Constraints. 1., calculate the inverse cell matrix U
    U = np.transpose( inv(L) )

    return L, U

#determine the volume of a cell given lengths and angles (a, b, c, alp, bet, gam)
def volume(length,angle):
    vlm3 = 1.
    for x in length:
        vlm3 *= x
    vlm1= 1.
    for x in angle:
        vlm1 -= np.cos(x)**2
    vlm2= 2.
    for x in angle:
        vlm2 *= np.cos(x)
    volume=vlm3*(vlm1+vlm2)**0.5
    return volume

sphcyl = False
lstar  = 0.0

# Boolean used to determine whether the unit cell under consideration is cubic or
# a general parallelipiped.
non_cubic = True

rotations = False
# Set to true if the rotational coordinates are to be given as unit vectors
unit_vector = False

lowestfile = open(file='./tet.dat', mode = 'r')

if(rotations):
    npart = int(int(lowestfile.readline())/2)
else:
    npart = int(int(lowestfile.readline()))

pos = np.zeros((npart,3))
if(rotations):
    ortn = np.zeros((npart,3))

skipline = lowestfile.readline().split()
for part in range(npart):
    pos[part] = [float(x) for x in lowestfile.readline().split()]

# print(pos)
if(rotations):
    for part in range(npart):
        q = [float(x) for x in lowestfile.readline().split()]
        if(unit_vector):
            ortn[part] = q[part]
        else:
            if(len(q)==4):
                ortn[part] = quat_to_unit(q)
            else:
                ortn[part] = aatouvec(q)

cell_lengths = [float(x) for x in lowestfile.readline().split()]
angles = to_radians([float(x) for x in lowestfile.readline().split()])
cellparam=np.concatenate([cell_lengths,angles],axis=0)
cm, icm = unit_cell(cellparam)

npartt = 4.0*npart

#convert the lines array into something useful
vlm = volume(cell_lengths,angles)
vlm_sphere = (4./3.)*math.pi*0.125
phi = (npartt*vlm_sphere) / vlm
print(phi)

gint = 999
grcut = 6.0
dgr = grcut / float(gint)

# grcnt = np.zeros((2,math.trunc(math.ceil(grcut/dgr)+1)))
grcnt = np.zeros(gint)
# print(grcnt)
grcutsq = grcut**2

na = [-10,-9,-8,-7,-6,-5,-4-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10]
nb = [-10,-9,-8,-7,-6,-5,-4-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10]
nc = [-10,-9,-8,-7,-6,-5,-4-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10]

# npartt = (npartt / 4.0)*len(na)*len(nb)*len(nc)

for a in na:
    for b in nb:
        for c in nc:
            n = (a,b,c)
            lat = np.matmul(n, cm)
            for i in range(npart):
                ri = pos[i]
                for j in range(npart):
                    if(a==0 and b==0 and c==0 and i==j):
                        continue                    
                    rj = pos[j]
                    rij = ri - rj + lat
                    dissq = np.matmul(rij,rij)
                    if(dissq < grcutsq):
                        gdis = math.sqrt(dissq)
                        gr = int( np.ceil(gdis/dgr) )
                        # if(a==0 and b==0 and c==0):
                        grcnt[gr] += 2.0
                        # else:
                            # grcnt[gr] += 2.0

g = []
y = []
for i in range(1,gint+1):
    g.append(dgr*(float(i)+0.5))
    da = ( (float(i+1)**3) - (float(i)**3) )*(dgr**3)
    da = da*math.pi*(4./3.)*phi
    y.append( grcnt[i-1]/(da*npart) )
    print(dgr*(float(i)+0.5), grcnt[i-1]/(da*npart))


ax.xaxis.set_tick_params(width=3.5, length=5.0, pad=10)
ax.yaxis.set_tick_params(width=3.5, length=5.0)
ax.xaxis.set_tick_params(width=2.5, length=3.5, which='minor',  pad=10)
ax.yaxis.set_tick_params(width=2.5, length=3.5, which='minor')

ax.xaxis.set_ticks(np.arange(0.0, grcut+0.1, 1.0))
# ax.yaxis.set_ticks(np.arange(0.0, 28.0, 2.0))
minor_ticks = np.arange(0.0, grcut+0.1, 0.5)
# minor_ticksy = np.arange(0.0, 28.0, 1.0)
ax.set_xticks(minor_ticks, minor=True)
# ax.set_yticks(minor_ticksy, minor=True)

ax.plot(g, y, 'o-', markersize=0, linewidth=1.0, color='black')
plt.show()