import math
import numpy as np
#from vedo import *
import random
import sys
import argparse
import copy
from pyknotid.spacecurves import Knot
from pyknotid.spacecurves import Link
import top_funcs

visitedList = []
ring_links  = []

def compute_linking_number(K,L):
    K = np.asarray(K)
    L = np.asarray(L)
    m = len(K)
    n = len(L)

    # lnm = np.zeros(m)
    # lnn = np.zeros(n)
    ln  = 0.

    for i in range(m):
        for j in range(n):
            a = L[j]       - K[i]
            b = L[j]       - K[(i+1)%m]
            c = L[(j+1)%n] - K[(i+1)%m]
            d = L[(j+1)%n] - K[i]

            an = np.linalg.norm(a)
            bn = np.linalg.norm(b)
            cn = np.linalg.norm(c)
            dn = np.linalg.norm(d)
            
            N1 = np.dot(a,np.cross(b,c))
            D1 = an*bn*cn + np.dot(a,b)*cn + np.dot(c,a)*bn + np.dot(b,c)*an

            N2 = np.dot(c,np.cross(d,a))
            D2 = cn*dn*an + np.dot(c,d)*an + np.dot(a,c)*dn + np.dot(d,a)*cn

            lni = ( np.arctan2(N1,D1) + np.arctan2(N2,D2) )
            ln  = ln + lni
            # lnm[i] = lnm[i] + lni
            # lnn[j] = lnn[j] + lni
    
    return ln / (2.*np.pi)#, lnm / (2.*np.pi), lnn / (2.*np.pi)

def compute_writhe(R):
    R = np.asarray(R)
    n = len(R)
    wr    = 0.#np.zeros(n)
    for i in range(n):
        for j in range(n):
            if(i==j or i==(j+1)%n or (i+1)%n==j or (i+1)%n==(j+1)%n):
                continue
            r1 = R[i]       - R[(j+1)%n]
            r2 = R[i]       - R[j]
            r3 = R[(i+1)%n] - R[(j+1)%n]
            r4 = R[(i+1)%n] - R[j]

    # ------------------------------------------------------------------------------
    #   Method from "Computation of Writhe in Modelling of Supercoiled DNA"
    # ------------------------------------------------------------------------------
            n1 = np.cross(r2,r1)
            n2 = np.cross(r1,r3)
            n3 = np.cross(r3,r4)
            n4 = np.cross(r4,r2)

            n1 = n1 / np.linalg.norm(n1)
            n2 = n2 / np.linalg.norm(n2)
            n3 = n3 / np.linalg.norm(n3)
            n4 = n4 / np.linalg.norm(n4)

            A = np.arcsin( np.dot(n1,n2) )
            B = np.arcsin( np.dot(n2,n3) )
            C = np.arcsin( np.dot(n3,n4) )
            D = np.arcsin( np.dot(n4,n1) )

            r12 = R[i] - R[(i+1)%n]
            r34 = R[j] - R[(j+1)%n]

            S = np.sign(np.dot(np.cross(r12,r34),r2))

            wr = wr + S*np.abs(A+B+C+D)
            # wr[j] = wr[j] + S*np.abs(A+B+C+D)
    # ------------------------------------------------------------------------------

    return wr/(4.*np.pi)

def depthFirst(graph, currentVertex, visited, nmin, nmax):
#   Check if any neighbours of the current vertex have previously been added to the chain.
#   If any have, and it isn't the initial node or the previous node it means there is a shorter path, 
#   and this chain should be ignored.
    for vertex in graph[currentVertex]:
        if vertex in visited:
            check1 = len(visited)==2 and vertex != visited[len(visited)-1]
            check2 = vertex != visited[0] and vertex != visited[len(visited)-1]
            if(check1 or check2):
                return
#   Add the current node to the chain of visited nodes
    visited.append(currentVertex)
#   If the length of the chain is greater than the max cycle length, end the recursive loop
    if(len(visited)>nmax):
        return
#   Check whether the initial node in the chain is a neighbour of the current node, if it is then
#   end the rescursion (by skipping this loop) as any potential rings that may form by adding more
#   nodes will not be primitive rings. 
    if(len(visited)<=2 or visited[0] not in graph[currentVertex] ):
        for vertex in graph[currentVertex]:
            if vertex not in visited:
                depthFirst(graph, vertex, visited.copy(), nmin, nmax)
#   If the length of the resulting chain is longer than the minimum cycle length, add it to 
#   the list of candidate cycles.
    if(len(visited)>=nmin):
#       Check if the chain is a cycle, if it is add it to the list of cycles.
        if( visited[0] in graph[visited[-1]] ):
            if( not visited[0] in graph[visited[-2]] ):
                global visitedList
                visitedList.append(visited)

#   Read conflink files as laid out by Professor Francesco Sciortino
def read_file(f):
    l = f.readline()
#   Read number of particles in the system
    npart = int(l.split()[0].split(',')[0])
#   Read in the box length
    boxl  = float(l.split()[1])
#   Initialise the variables to be populated by reading the conflink file.
    r  = np.zeros((npart,3))
    nb = np.zeros(npart)
    g  = {}
#   Read in the positions of the particles
    for i in range(npart):
        l = f.readline().split()
        r[i] = float(l[0].split(',')[0]), float(l[1].split(',')[0]), float(l[2])
#   Populate the dictionary defining the graph network of the system under consideration.
    for i in range(npart):
        l     = f.readline().split()
        nb[i] = int(l[1])
        l_nbs = f.readline().split()
        nbs   = []
        for j in l_nbs:
            nbs.append(int(j.split(',')[0])-1)
        g[i] = nbs

    return npart, boxl, r, nb, g

def main(args):

    ntraj       = args.ntraj #100
    nprint      = args.nprint #10
    ring_min    = args.min #4
    ring_max    = args.max #16
    baset       = args.base
    watert      = args.water

    if(watert):
        dcut = 0.4**2
    else:
        dcut = 9.

    path2file   = './conflink.dat'
    f           = open(path2file,mode='r')
    count_file  = open('counts.out',mode='w+')
    ringid_file = open('ring_id.out',mode='w+')
    wrn_file    = open('wr.dat',mode='w+')
    wrnd_file   = open('wr_dist.dat',mode='w+')
    unknot_file = open('wr_unknot.dat',mode='w+')
    #hc_file     = open('helicity.xyz',mode='w+')
    printt      = False
# ===================================================================================
# ===================================================================================
#   BEGIN IDENTIFICATION OF RINGS IN THE NETWORK
# ===================================================================================
# ===================================================================================
    for ni in range(ntraj):
        global visitedList
        visitedList = []
        npart, boxl, r, nb, graph = read_file(f)
        # copy_graph = copy.deepcopy(graph)
    #   Loop over all particles in order to extract the chains of particles that correspond to rings.
        for i in range(npart):
        #   Call DFS method to extract all chains of particles that correspond to rings, where particle i
        #   is the first particle in the chain. 
            depthFirst(graph, i, [], ring_min, ring_max)
        #   Remove all connections to particle i so that rings involving particle i are not counted again.
        #   Note however that all rings will be counted twice as there are two possible ways of goind round a
        #   ring (i.e., forwards or backwards).
            del graph[i]
            for l in graph.values():
                if i in l: l.remove(i)

        counts = {}
        prunedList = []
        for iv in range(len(visitedList)-1):
            vring1 = set(visitedList[iv])
            for jv in range(iv+1,len(visitedList)):
                vring2 = set(visitedList[jv])
                if(vring1==vring2):
                    prunedList.append(visitedList[iv])
        visitedList = prunedList

        visitedList = list(visitedList)
        visitedList.sort(key=len)

#   -------------------------------------------------------------------------------------
#       Remove percolating chains from the list of rings
#   -------------------------------------------------------------------------------------
        rings_pos = []
        ring_list = []
        ring_com  = []
        for _,v in enumerate(visitedList):
            rp = []

            for i in v:
                rp.append(r[i])
            max_d, rp, com = top_funcs.perc_check(rp,boxl,dcut)

            if(max_d <= dcut):
                rings_pos.append(rp)
                ring_list.append(v)
                ring_com.append(com)

        visitedList = ring_list
#   -------------------------------------------------------------------------------------

#   -----------------------------------------------------------------------------------------
#       If considering only shortest path rings, remove rings that are a linear combination
#       of pairs and triplets of smaller rings.
#   -----------------------------------------------------------------------------------------
        if(baset):
            base_rings = copy.deepcopy(visitedList)
            base_pos   = copy.deepcopy(rings_pos)
            base_com   = copy.deepcopy(ring_com)
            
            nr = len(visitedList)
        
        #   Prune list of rings to keep only those that are base loops.
            for rip in range(nr-1):
                ring_i = visitedList[rip]
                
                for rjp in range(rip+1,nr):
                    ring_j = visitedList[rjp]
                    inter = list(set.intersection(set(ring_i), set(ring_j)))
                    
                    if(len(inter)>0):
                        nodes_i = copy.deepcopy(ring_i)
                        nodes_j = copy.deepcopy(ring_j)
                        for node in inter:
                            nodes_i.remove(node)
                            nodes_j.remove(node)

                        lb = len(base_rings)-1
                        if(len(inter)>2):
                            for ki,ring_k in enumerate(reversed(base_rings)):
                                if(ring_k != ring_i and ring_k != ring_j):
                                    if(len(ring_k)>len(ring_i) and len(ring_k)>len(ring_j)):
                                        if(len(set(ring_k) & set(inter)) >= 2):
                                            if( all(x in ring_k for x in nodes_i) ):
                                                if( all(x in ring_k for x in nodes_j) ):
                                                    remove_index = lb-ki
                                                    del base_rings[remove_index]
                                                    del base_pos[remove_index]
                                                    del base_com[remove_index]
                                                    break
                        else:
                            test_ring = nodes_i + nodes_j
                            for ki,ring_k in enumerate(reversed(base_rings)):
                                if(len(set(test_ring) & set(ring_k)) == len(test_ring)):
                                    remove_index = lb-ki
                                    del base_rings[remove_index]
                                    del base_pos[remove_index]
                                    del base_com[remove_index]
                                    break
                            
                        for rkp in range(rjp+1,nr):
                            ring_k = visitedList[rkp]
                            inter_ik = list(set.intersection(set(ring_i), set(ring_k)))
                            inter_jk = list(set.intersection(set(ring_j), set(ring_k)))
                            
                            if( len(inter_ik)>0 and len(inter_jk)>0 ):
                                inter_tot = list(set(inter + inter_ik + inter_jk))
                                nodes_i = copy.deepcopy(ring_i)
                                nodes_j = copy.deepcopy(ring_j)
                                nodes_k = copy.deepcopy(ring_k)
                                
                                for node in inter_tot:
                                    if(node in nodes_i):
                                        nodes_i.remove(node)
                                    if(node in nodes_j):
                                        nodes_j.remove(node)
                                    if(node in nodes_k):
                                        nodes_k.remove(node)
                                
                                lb = len(base_rings)-1
                                for li,ring_l in enumerate(reversed(base_rings)):
                                    if(ring_l != ring_i and ring_l != ring_j and ring_l != ring_k):
                                        if(len(ring_l)>len(ring_i) and len(ring_l)>len(ring_j) and len(ring_l)>len(ring_k)):
                                            if(len(set(ring_l) & set(inter_tot)) >= 3):
                                                if( len(set(ring_l) & set(nodes_i)) >= len(nodes_i)-1 ):
                                                    if( len(set(ring_l) & set(nodes_j)) >= len(nodes_j)-1 ):
                                                        if( len(set(ring_l) & set(nodes_k)) >= len(nodes_k)-1 ):
                                                            remove_index = lb-li
                                                            del base_rings[remove_index]
                                                            del base_pos[remove_index]
                                                            del base_com[remove_index]
                                                
            visitedList = copy.deepcopy(base_rings)
            rings_pos   = copy.deepcopy(base_pos)
            ring_com    = copy.deepcopy(base_com)
#   -----------------------------------------------------------------------------------------

#   -----------------------------------------------------------------------------------------

    #   Print data related to ring statistics and open files related to knots    
        for i in range(ring_min, ring_max+1):
            counts[i] = 0
        for v in visitedList:
            counts[len(v)] += 1
        for i in range(ring_min, ring_max+1):
            count_file.write(str(counts[i]) + ' ')
        count_file.write('\n')
        count_file.flush()

        # wr1_tot     = 0.0
        # wr2_tot     = 0.0
        # wr3_tot     = 0.0
        unknot_sum  = 0
        knot_sum    = 0
        wr_tot      = 0.0
        gln_sum     = 0.0

        if(printt and (ni==0 or (ni+1)%nprint==0)):
            ring_file = open('rings'+str(ni+1)+'.xyz',mode='w+')
            link_file = open('links'+str(ni+1)+'.xyz',mode='w+')
            wr1_file  = open('trefoil'+str(ni+1)+'.xyz',mode='w+')
            wr2_file  = open('figeigt'+str(ni+1)+'.xyz',mode='w+')
            # for vi,(v,vp) in enumerate(zip(visitedList,rings_pos)):

                # ring_file.write(str(npart)+'\n')
                # ring_file.write(str(len(v)) + ' ' + str(v)+'\n')
                
                # for i in range(npart):
                #     if(i not in v):
                #         ri = r[i]
                #         ring_file.write('C ' + str(ri[0]) + ' ' + str(ri[1]) + ' ' + str(ri[2]) + '\n')
                
                # for ri in vp:
                #     ring_file.write('O ' + str(ri[0]) + ' ' + str(ri[1]) + ' ' + str(ri[2]) + '\n')

            # ring_file.close()
        if(printt):
            wr3_file = open('bigknot'+str(ni+1)+'.xyz',mode='w+')

# ===================================================================================
# ===================================================================================
#   BEGIN IDENTIFICATION OF RINGS IN THE NETWORK
# ===================================================================================
# ===================================================================================
        
    # -----------------------------------------------------------
    #   Identify individual rings that are knotted
    # -----------------------------------------------------------

        knot_rings = copy.deepcopy(visitedList)
        knot_pos   = copy.deepcopy(rings_pos)
        knot_com   = copy.deepcopy(ring_com)

        lb = len(rings_pos)-1

        # hc_i  = np.zeros(npart)
        
        for rip,rpi in enumerate(reversed(rings_pos)):
            
            remove_index = lb-rip
            ring_i       = visitedList[remove_index]
            
            try:
                k1  = Knot(np.asarray(rpi),verbose=False)
                gc1 = k1.gauss_code()
                gc1.simplify(one=True,two=True,one_extended=True)
            except IndexError:
                continue

            wr = top_funcs.calc_writhe(np.asarray(rpi))
            wr_tot = wr_tot + abs(wr)

            if(printt and (ni==0 or (ni+1)%nprint==0)):
                ring_file.write(str(len(rpi))+'\n')
                ring_file.write(str(wr) + ', ' + str(ring_i)+'\n')
                
                for ri in rpi:
                    ring_file.write('O ' + str(ri[0]) + ' ' + str(ri[1]) + ' ' + str(ri[2]) + '\n')

            if(str(gc1)!='----'):
                
                knot_id = k1.identify()

                if(len(knot_id)>0):
                    strk    = str(knot_id[0])
                else:
                    strk = '<Knot 0_1>'
                
                if(strk=='<Knot 0_1>'):
                    unknot_sum = unknot_sum + 1
                    unknot_file.write(str(len(rpi))+' '+str(abs(wr))+'\n')
                    unknot_file.flush()

                else:#if(strk!='<Knot 0_1>'):
                    
                    knot_sum = knot_sum + 1
                    # wr = top_funcs.calc_writhe(np.asarray(rpi))
                    # wr_tot = wr_tot + abs(wr)
                    wrnd_file.write(strk[6:9]+' '+str(abs(wr))+'\n')
                    wrnd_file.flush()

                #   Remove rings from the list to prevent overcounting of knots
                    del knot_rings[remove_index] 
                    del knot_pos[remove_index]
                    del knot_com[remove_index]

                    if(printt and (ni==0 or (ni+1)%nprint==0)):
                        if(strk=='<Knot 3_1>'):
                            wr1_file.write(str(len(ring_i))+'\n')
                            wr1_file.write(str(wr)+ '; ' + str(knot_id)+ '; ' + str(gc1)+ '; ' + str(ring_i) + '\n')
                            for ri in rpi:
                                wr1_file.write('O ' + str(ri[0]) + ' ' + str(ri[1]) + ' ' + str(ri[2]) + '\n')
                            # wr1_file.flush()
                        elif(strk=='<Knot 4_1>'):
                            wr2_file.write(str(len(ring_i))+'\n')
                            wr2_file.write(str(wr)+ '; ' + str(knot_id)+ '; ' + str(gc1)+ '; ' + str(ring_i) + '\n')
                            for ri in rpi:
                                wr2_file.write('O ' + str(ri[0]) + ' ' + str(ri[1]) + ' ' + str(ri[2]) + '\n')
                            # wr2_file.flush()
                        else:
                            wr3_file.write(str(len(ring_i))+'\n')
                            wr3_file.write(str(wr)+ '; ' + str(knot_id)+ '; ' + str(gc1)+ '; ' + str(ring_i) + '\n')
                            for ri in rpi:
                                wr3_file.write('O ' + str(ri[0]) + ' ' + str(ri[1]) + ' ' + str(ri[2]) + '\n')
                            # wr3_file.flush()
                    else:
                        if(printt and strk!='<Knot 3_1>' and strk!='<Knot 4_1>'):
                            wr3_file.write(str(len(ring_i))+'\n')
                            wr3_file.write(str(wr)+ '; ' + str(knot_id)+ '; ' + str(gc1)+ '; ' + str(ring_i) + '\n')
                            for ri in rpi:
                                wr3_file.write('O ' + str(ri[0]) + ' ' + str(ri[1]) + ' ' + str(ri[2]) + '\n')
                            # wr3_file.flush()
            else:
                unknot_sum = unknot_sum + 1
                unknot_file.write(str(len(rpi))+' '+str(abs(wr))+'\n')
                unknot_file.flush()
        
        if(printt and (ni==0 or (ni+1)%nprint==0)):
            ring_file.close()
    
    # -----------------------------------------------------------------------
    #  Identify pairs of joint rings that are knotted, such as theta-curves
    # -----------------------------------------------------------------------
        visitedList = copy.deepcopy(knot_rings)
        rings_pos   = copy.deepcopy(knot_pos)
        ring_com    = copy.deepcopy(knot_com)
        
        for rip in range(len(rings_pos)-1):
            ring_i = visitedList[rip]
            rpi    = rings_pos[rip]
            comi   = ring_com[rip]

            for rjp in range(rip+1,len(rings_pos)):
                
                ring_j = visitedList[rjp]
                inter = list(set.intersection(set(ring_i), set(ring_j)))
                
                if(len(inter)>1):
                    rpj    = rings_pos[rjp]
                    comj   = ring_com[rjp]
                    
                # -----------------------------------------------------------
                #   Use minimum image convention to identify nearest image of
                #   particles in Ring J to Ring I, using the COM of Ring I.
                # -----------------------------------------------------------
                    rpj = top_funcs.ring_image(comi,rpj,comj,boxl)

                    pint_i = ring_i.index(inter[0])
                    pint_j = ring_j.index(inter[0])

                    ring_ij = np.asarray(rpi[pint_i]) - np.asarray(rpj[pint_j])
                    max_d   = np.dot(ring_ij,ring_ij)
                # -----------------------------------------------------------

                # -----------------------------------------------------------
                #   Check that shared particles are next to one another in 
                #   each of the rings, otherwise the minimum image convention
                #   has failed.
                # -----------------------------------------------------------
                    if(max_d <= 1.):
                        rp1   = []
                        rp2   = []

                    # -----------------------------------------------------------
                    #   Determine proper sequence of particles to compute the 
                    #   Writhe and define a knot.
                    # -----------------------------------------------------------
                        anchor = 0
                        for int_i in inter:
                            mn1 = (ring_i.index(int_i) - 1)%len(ring_i)
                            if(ring_i[mn1] not in inter):
                                anchor = int_i
                                break
                        
                        new_inter = []
                        for int_i in range(len(inter)):
                            inter_index = (ring_i.index(anchor)+int_i)%len(ring_i)
                            new_inter.append( ring_i[inter_index] )
                        
                        if(len(set.intersection(set(inter), set(new_inter)))==len(inter)):
                            inter = copy.deepcopy(new_inter)

                            first_node = ring_i.index(inter[0])
                            reverse1 = False
                            if(ring_i[(first_node+1)%len(ring_i)] in inter):
                                reverse1 = True
                            ring1 = []
                            for i in range(1,len(ring_i)-len(inter)+1):
                                if(reverse1):
                                    index = (first_node-i)%len(ring_i)
                                else:
                                    index = (first_node+i)%len(ring_i)
                                ring1.append(ring_i[index])
                                rp1.append(rpi[index])

                            first_node = ring_j.index(inter[-1])
                            reverse2 = False
                            if(ring_j[(first_node+1)%len(ring_j)] in inter):
                                reverse2 = True
                            ring2 = []
                            for i in range(1,len(ring_j)-len(inter)+1):
                                if(reverse2):
                                    index = (first_node-i)%len(ring_j)
                                else:
                                    index = (first_node+i)%len(ring_j)
                                ring2.append(ring_j[index])
                                rp2.append(rpj[index])

                            path = []
                            for path_ind in inter:
                                p1   = ring_i.index(path_ind)
                                path.append(rpi[p1])
                    # -----------------------------------------------------------

                            if(inter[0] in ring1 or inter[1] in ring1 or inter[0] in ring2 or inter[1] in ring2):
                                continue
                            else:
                            # -----------------------------------------------------------
                            #   Define potential knot and compute simplified Gauss code
                            # -----------------------------------------------------------
                                id_path   = [inter[0]]+ring1+[inter[-1]]+ring2
                                ring_path = [path[0]]+rp1+[path[-1]]+rp2
                                try:
                                    k  = Knot(np.asarray(ring_path),verbose=False)
                                    gc = k.gauss_code()
                                    gc.simplify(one=True,two=True,one_extended=True)
                                except IndexError:
                                    continue
                            # -----------------------------------------------------------

                            #   If the ring is not an Unknot, continue to perform further analysis
                                if(str(gc)!='----'):
                                # -----------------------------------------------------------------------------
                                #  Again, if the ring is not an Unknot, continue to perform further analysis.
                                #  Also, if Ring I has already been identified as part of a knot of the same
                                #  type do not count that knot again.
                                # -----------------------------------------------------------------------------
                                    knot_id = k.identify()

                                    if(len(knot_id)>0):
                                        strk    = str(knot_id[0])
                                    else:
                                        strk = '<Knot 0_1>'

                                    if(strk=='<Knot 0_1>'):
                                        unknot_sum = unknot_sum + 1
                                    else:
                                    # -----------------------------------------------------------
                                    #   Compute the Writhe for the Knot and print data to files
                                    # -----------------------------------------------------------
                                        knot_sum = knot_sum + 1
                                        wr = top_funcs.calc_writhe(np.asarray(ring_path))
                                        wr_tot = wr_tot + abs(wr)
                                        
                                        wrnd_file.write(strk[6:9]+' '+str(abs(wr))+'\n')
                                        wrnd_file.flush()

                                        if(printt and (ni==0 or (ni+1)%nprint==0)):
                                            if(strk=='<Knot 3_1>'):
                                                wr1_file.write(str(len(rp1)+len(rp2)+len(path))+'\n')
                                                wr1_file.write(str(wr)+ '; ' + strk + '; ' + str(gc)+ '; ' + str(ring1) + '; ' + str(ring2) + '; ' + str(inter) + '\n')
                                                for ri in rp1:
                                                    wr1_file.write('O ' + str(ri[0]) + ' ' + str(ri[1]) + ' ' + str(ri[2]) + '\n')
                                                for rj in rp2:
                                                    wr1_file.write('C ' + str(rj[0]) + ' ' + str(rj[1]) + ' ' + str(rj[2]) + '\n')
                                                for rl in path:
                                                    wr1_file.write('N ' + str(rl[0]) + ' ' + str(rl[1]) + ' ' + str(rl[2]) + '\n')
                                                # wr1_file.flush()
                                            elif(strk=='<Knot 4_1>'):
                                                wr2_file.write(str(len(rp1)+len(rp2)+len(path))+'\n')
                                                wr2_file.write(str(wr)+ '; ' + strk + '; ' + str(gc)+ '; ' + str(ring1) + '; ' + str(ring2) + '; ' + str(inter) + '\n')
                                                for ri in rp1:
                                                    wr2_file.write('O ' + str(ri[0]) + ' ' + str(ri[1]) + ' ' + str(ri[2]) + '\n')
                                                for rj in rp2:
                                                    wr2_file.write('C ' + str(rj[0]) + ' ' + str(rj[1]) + ' ' + str(rj[2]) + '\n')
                                                for rl in path:
                                                    wr2_file.write('N ' + str(rl[0]) + ' ' + str(rl[1]) + ' ' + str(rl[2]) + '\n')
                                                # wr2_file.flush()
                                            else:
                                                wr3_file.write(str(len(rp1)+len(rp2)+len(path))+'\n')
                                                wr3_file.write(str(wr)+ '; ' + strk + '; ' + str(gc)+ '; ' + str(ring1) + '; ' + str(ring2) + '; ' + str(inter) + '\n')
                                                for ri in rp1:
                                                    wr3_file.write('O ' + str(ri[0]) + ' ' + str(ri[1]) + ' ' + str(ri[2]) + '\n')
                                                for rj in rp2:
                                                    wr3_file.write('C ' + str(rj[0]) + ' ' + str(rj[1]) + ' ' + str(rj[2]) + '\n')
                                                for rl in path:
                                                    wr3_file.write('N ' + str(rl[0]) + ' ' + str(rl[1]) + ' ' + str(rl[2]) + '\n')
                                                # wr3_file.flush()
                                        else:
                                            if(printt and strk!='<Knot 3_1>' and strk!='<Knot 4_1>'):
                                                wr3_file.write(str(len(rp1)+len(rp2)+len(path))+'\n')
                                                wr3_file.write(str(wr)+ '; ' + strk + '; ' + str(gc)+ '; ' + str(ring1) + '; ' + str(ring2) + '; ' + str(inter) + '\n')
                                                for ri in rp1:
                                                    wr3_file.write('O ' + str(ri[0]) + ' ' + str(ri[1]) + ' ' + str(ri[2]) + '\n')
                                                for rj in rp2:
                                                    wr3_file.write('C ' + str(rj[0]) + ' ' + str(rj[1]) + ' ' + str(rj[2]) + '\n')
                                                for rl in path:
                                                    wr3_file.write('N ' + str(rl[0]) + ' ' + str(rl[1]) + ' ' + str(rl[2]) + '\n')
                                                # wr3_file.flush()
                                else:
                                    unknot_sum = unknot_sum + 1
                elif(len(inter)==0):
                # ---------------------------------------------------------------------------
                # Identify if the two disjoint rings are linked.
                # ---------------------------------------------------------------------------   
                    rpj    = rings_pos[rjp]
                    comj   = ring_com[rjp]
                    
                    rpj = top_funcs.ring_image(comi,rpj,comj,boxl)

                    lk = top_funcs.calc_linking_number(np.asarray(rpi),np.asarray(rpj))
                    lk = np.abs(lk)

                    if(lk>1.0):
                        gln_sum = gln_sum + int(lk/2)
                            
                        if(printt and (ni==0 or (ni+1)%nprint==0)):
                            link_file.write(str(len(rpi)+len(rpj))+'\n')
                            link_file.write(str(lk) + '; ' + str(ring_i) + '; ' + str(ring_j) + '\n')
                            for ri in rpi:
                                link_file.write('O ' + str(ri[0]) + ' ' + str(ri[1]) + ' ' + str(ri[2]) + '\n')
                            for rj in rpj:
                                link_file.write('C ' + str(rj[0]) + ' ' + str(rj[1]) + ' ' + str(rj[2]) + '\n')
                            # link_file.flush()
                    continue
                
                else:
                    continue

        n_rings = sum(counts.values())
        wrn_file.write(str(gln_sum)+' '+str(wr_tot)+' '+str(n_rings)+' '+str(gln_sum/n_rings)+' '+str(wr_tot/n_rings)+' '+str((gln_sum+wr_tot)/n_rings)+' '+str(knot_sum)+' '+str(unknot_sum)+'\n')
        wrn_file.flush()
        
        if(printt and (ni==0 or (ni+1)%nprint==0)):
            wr1_file.close()
            wr2_file.close()
            wr3_file.close()
            link_file.close()

        # hc_file.write(str(npart)+'\n')
        # hc_file.write(str(wr_tot)+'\n')
        # for npi in range(npart):
        #     ri = r[npi]
        #     hc_file.write('C '+str(ri[0])+' '+str(ri[1])+' '+str(ri[2])+' '+str(hc_i[npi]))
    unknot_file.close()
# =============================================================
# Define parameters which can be read in from the command line
# =============================================================
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
# Paramters defining the diamond system
    parser.add_argument('-ntraj',    '--ntraj',    type=int,            default=1,     help='integer giving the number of simulation snapshots')
    parser.add_argument('-nprint',   '--nprint',   type=int,            default=1,     help='integer giving the frequency to print visualisations of the rings.')
    parser.add_argument('-min',      '--min',      type=int,            default=4,     help='integer giving the minimum ring size to extract.')
    parser.add_argument('-max',      '--max',      type=int,            default=9,     help='integer giving the maximum ring size to extract.')
    parser.add_argument('-base',     '--base',     action='store_true', default=False, help='Boolean set to true if using rings in base set only.')
    parser.add_argument('-water',    '--water',    action='store_true', default=False, help='Boolean set to true if analysing water configurations.')
    args = parser.parse_args()
    main(args)