# https://graphics.stanford.edu/papers/fastlinkingnumbers/assets/QuJames2021FastLinkingNumbers_Optimized400dpi.pdf
# https://arxiv.org/pdf/1912.13121.pdf
# https://www.jstage.jst.go.jp/article/nolta/4/1/4_104/_pdf/-char/en

import numpy as np
import sys
import argparse
import copy

visitedList = []
ring_links  = []

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
    watertt     = args.water

    path2file             = './conflink.dat'
    f                     = open(path2file,mode='r')
    count_file            = open('counts.out',mode='w+')
    ringid_file           = open('ring_id.out',mode='w+')

    for ni in range(ntraj):
        global visitedList
        visitedList = []
        npart, boxl, r, nb, graph = read_file(f)
        copy_graph = copy.deepcopy(graph)
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
        
        ringid_file.write("Snapshot: " + str(ni) + ", " + str(len(visitedList)) + " unique rings found." + '\n')
        for v in visitedList:
            ringid_file.write(str(len(v)) + ' ' + str(v)+'\n')

        rings_pos = []
        ring_list = []
        ring_com  = []
        for vi,v in enumerate(visitedList):
            rp = []
            r0 = r[v[0]]
            
            for i in v:
                rij = r0 - r[i]
                rij = rij - (np.rint(rij / boxl) * boxl)
                ri  = r0 - rij
                rp.append([ri[0],ri[1],ri[2]])

            max_d = 0.
            for i in range(len(rp)):
                ring_ij = np.asarray(rp[i]) - np.asarray(rp[(i+1)%(len(rp))])
                ring_d  = np.dot(ring_ij,ring_ij)
                if(ring_d > max_d):
                    max_d = ring_d

            if(watertt):
                if(max_d <= 0.4**2):
                    rings_pos.append(rp)
                    ring_list.append(v)
            else:
                if(max_d <= 9.):
                    rings_pos.append(rp)
                    ring_list.append(v)
        
        visitedList = ring_list

        for rip in range(len(rings_pos)):
            rpi    = rings_pos[rip]
            comi   = np.zeros(3)
            for ri in rpi:
                comi = comi + ri
            comi = comi / len(rpi)
            ring_com.append(comi)

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
        
        for i in range(ring_min, ring_max+1):
            counts[i] = 0
        for v in visitedList:
            counts[len(v)] += 1
        for i in range(ring_min, ring_max+1):
            count_file.write(str(counts[i]) + ' ')
        count_file.write('\n')
        count_file.flush()

        if(ni==0 or (ni+1)%nprint==0):
            ring_file = open('rings'+str(ni+1)+'.xyz',mode='w+')
            for vi,(v,vp) in enumerate(zip(visitedList,rings_pos)):

                ring_file.write(str(npart)+'\n')
                ring_file.write(str(len(v)) + ' ' + str(v)+'\n')
                
                for i in range(npart):
                    if(i not in v):
                        ri = r[i]
                        ring_file.write('C ' + str(round(ri[0],4)) + ' ' + str(round(ri[1],4)) + ' ' + str(round(ri[2],4)) + '\n')
                
                for ri in vp:
                    ring_file.write('O ' + str(round(ri[0],4)) + ' ' + str(round(ri[1],4)) + ' ' + str(round(ri[2],4)) + '\n')

            ring_file.close()

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
    parser.add_argument('-water',    '--water',    action='store_true', default=False, help='Boolean set to true if analysing water configurations.')
    args = parser.parse_args()
    main(args)