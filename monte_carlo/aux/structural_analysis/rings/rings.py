import numpy as np
import sys
import argparse

visitedList = []

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

    path2file   = './conflink.dat'
    f           = open(path2file,mode='r')
    count_file  = open('counts.out',mode='w+')
    ringid_file = open('ring_id.out',mode='w+')

    for ni in range(ntraj):
        global visitedList
        visitedList = []
        npart, boxl, r, nb, graph = read_file(f)
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
        visitedList = {tuple(item) for item in map(sorted, visitedList)}
        for i in range(ring_min, ring_max+1):
            counts[i] = 0
        for v in visitedList:
            counts[len(v)] += 1
        for i in range(ring_min, ring_max+1):
            count_file.write(str(counts[i]) + ' ')
        count_file.write('\n')
        print(ni,counts)

        visitedList = list(visitedList)
        visitedList.sort(key=len)
        ringid_file.write("Snapshot: " + str(ni) + ", " + str(len(visitedList)) + " unique rings found." + '\n')
        for v in visitedList:
            ringid_file.write(str(len(v)) + ' ' + str(v)+'\n')

        if(ni==0 or (ni+1)%nprint==0):
            ring_file = open('rings'+str(ni+1)+'.xyz',mode='w+')
            for v in visitedList:

                ring_file.write(str(npart)+'\n')
                ring_file.write(str(len(v)) + ' ' + str(v)+'\n')
                
                for i in range(npart):
                    if(i not in v):
                        ri = r[i]
                        ring_file.write('C ' + str(round(ri[0],4)) + ' ' + str(round(ri[1],4)) + ' ' + str(round(ri[2],4)) + '\n')
                
                r0 = r[v[0]]
                for i in v:
                    rij = r0 - r[i]
                    rij = rij - (np.rint(rij / boxl) * boxl)
                    ri  = r0 - rij
                    ring_file.write('O ' + str(round(ri[0],4)) + ' ' + str(round(ri[1],4)) + ' ' + str(round(ri[2],4)) + '\n')

            ring_file.close()

# =============================================================
# Define parameters which can be read in from the command line
# =============================================================
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
# Paramters defining the diamond system
    parser.add_argument('-ntraj',  '--ntraj',  type=int, default=1, help='integer giving the number of simulation snapshots')
    parser.add_argument('-nprint', '--nprint', type=int, default=1, help='integer giving the frequency to print visualisations of the rings.')
    parser.add_argument('-min',    '--min',    type=int, default=4, help='integer giving the minimum ring size to extract.')
    parser.add_argument('-max',    '--max',    type=int, default=9, help='integer giving the maximum ring size to extract.')
    args = parser.parse_args()
    main(args)