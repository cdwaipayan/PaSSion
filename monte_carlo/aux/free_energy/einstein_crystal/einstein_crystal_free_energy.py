"""

"""
import argparse
import numpy as np
import vegas
import math

def rigid_ref(nsites):

    refsites = np.zeros( (nsites,3) )

    if(nsites == 1):
        refsites[:][0] = [0.0, 0.0, 1.0]

    elif(nsites == 2):
        refsites[:][0] = [0.0, 0.0, 1.0]
        refsites[:][1] = [0.0, 0.0,-1.0]

    elif(nsites == 3):
        refsites[:][0] = [ 0.0,          0.0,           1.0 ]
        refsites[:][1] = [ 0.0,  (math.sqrt(3.0)/2.0), -0.5 ]
        refsites[:][2] = [ 0.0, -(math.sqrt(3.0)/2.0), -0.5 ]

    elif(nsites == 4):
        refsites[:][0] = [ 1.0/math.sqrt(3.0),  1.0/math.sqrt(3.0),  1.0/math.sqrt(3.0)]
        refsites[:][1] = [-1.0/math.sqrt(3.0), -1.0/math.sqrt(3.0),  1.0/math.sqrt(3.0)]
        refsites[:][2] = [ 1.0/math.sqrt(3.0), -1.0/math.sqrt(3.0), -1.0/math.sqrt(3.0)]
        refsites[:][3] = [-1.0/math.sqrt(3.0),  1.0/math.sqrt(3.0), -1.0/math.sqrt(3.0)]

    elif(nsites == 6):
        refsites[:][0] = [ 1.0, 0.0, 0.0]
        refsites[:][1] = [ 0.0, 1.0, 0.0]
        refsites[:][2] = [ 0.0, 0.0, 1.0]
        refsites[:][3] = [-1.0, 0.0, 0.0]
        refsites[:][4] = [ 0.0,-1.0, 0.0]
        refsites[:][5] = [ 0.0, 0.0,-1.0]

    return refsites


def einstein_crystal_or(nsites, rbref, new_rbsites):
    
    energy = 0.0

    ja  = 0
    jb  = 0
    jb2 = 0

    for j1 in range(nsites):
#   calculate the angle subtended by the position of reference site 1 and each of the new sites. 
#   additionally, we make sure to account for rounding errors in double precision floating point numbers 
        thetaa = math.acos(min(max(np.matmul( rbref[0], new_rbsites[j1] ), -1.0), 1.0))
    #   keep track of which site gives the smallest angle
        if(j1 == 0 or thetaa < thetaa_min):
            ja = j1
            thetaa_min = thetaa

        if(nsites > 2):
#   calculate the angle subtended by the position of reference site 2 and each of the new sites. 
#   additionally, we make sure to account for rounding errors in double precision floating point numbers 
            thetab = math.acos(min(max(np.matmul( rbref[1],new_rbsites[j1] ), -1.0), 1.0))
        #   keep track of the smallest angle
            if(j1 == 0 or thetab < thetab_min):

                if(j1 >= 1):
                    jb2 = jb
                    thetab_min2 = thetab_min

                jb = j1
                thetab_min  = thetab

            elif(j1 == 1 or thetab < thetab_min2):
                jb2 = j1
                thetab_min2 = thetab
    
#   if the same rigid body site in the new configuration provides the smallest angle 
#   between both of the reference sites, reject the move.
    if(ja == jb):
        jb = jb2
        thetab_min = thetab_min2

    if(nsites<=2):
#   particles have uniaxial symmetry and so we only need to consider one of the sites.
        energy = math.sin(thetaa_min)**2
    else:
#   particles do not have uniaxial symmetry and so we need to consider at least two sites to constrain
#   the orientation of the particles.
        if(nsites==3): # Particles with D3h symmetry 
            energy = math.sin(3.0*thetaa_min/2.0)**2 + math.sin(thetab_min)**2 

        elif(nsites==4 or nsites==6): # Particles with Oh or Td symmetry
            energy = math.sin(thetaa_min)**2 + math.sin(thetab_min)**2 

    return energy


def angular_partition_function(x):

    rm = np.zeros( (3,3) )
    new_rbsites = np.zeros( (nsites, 3) )

    phi   = x[0] 
    theta = math.acos(x[1]) 
    psi   = x[2]

#   Calculate the corresponding rotation matrix using the Z-X-Z convention
    c1 = math.cos(phi); c2 = math.cos(theta); c3 = math.cos(psi)
    s1 = math.sin(phi); s2 = math.sin(theta); s3 = math.sin(psi)

    rm[0][:] = [ c1*c3-s1*c2*s3, -c1*s3-s1*c2*c3,  s1*s2]
    rm[1][:] = [ s1*c3+c1*c2*s3, -s1*s3+c1*c2*c3, -c1*s2]
    rm[2][:] = [    s2*s3,           s2*c3,         c2  ]

#   rotate the rigid body sites on the reference particle
    for j2 in range(nsites):
        new_rbsites[:][j2] = np.matmul(rm, rbref[:][j2])
        
#   calculate the energy associated with rotating the particles
    uor =  einstein_crystal_or(nsites, rbref, new_rbsites)
    
    return math.exp(-lamor*uor)

# =============================
# MAIN BODY OF SCRIPT
# =============================
pi    = math.pi
twopi = 2.0*math.pi

# Parameter to change
nitn  = 100
neval = 2500000

nsites = 4
lamor  = 10000

rbref  = rigid_ref(nsites) 

print( 'nitn = %s  neval = %s \n' % (nitn, neval) )
print( 'lamor = %s \n' % (lamor) )

integ = vegas.Integrator([[0, twopi], [-1, 1], [0, twopi]])
result = integ(angular_partition_function, nitn=nitn, neval=neval)

print(result.summary())
print('result = %s    Q = %.2f' % (result, result.Q))
print('Orientational Free Energy per Particle = %s' % -math.log(result.mean / (8.0*pi**2)) )
