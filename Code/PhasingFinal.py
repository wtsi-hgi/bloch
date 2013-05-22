# Copyright (c) 2012 Genome Research Ltd. Author: Michelle Parker <mp18@sanger.ac.uk>
# This file is part of BLOCH.
# BLOCH is free software: you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation;
# either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

#Import required python packages
import csv
import networkx as nx
import random
import numpy as np
import itertools
import math

class HMM:
    def __init__(self, T, gt):
        self.T = T
        self.gt = gt
        
    #Haploid initial state probabilities. allele: allele number.
    def hapinitial(self, allele):    
        #Edge counts are used. Count of edge/Total count for all edges.
        for i in self.T[0].out_edges(1, data=True):        
            if i[2]['a'] == allele:           
                return float(i[2]['weight'])/float(sum(self.T[0].node[1]['hap'].itervalues()))

    #Diploid initial state probabilities. a,b: allele number
    def dipinitial(self, a,b):    
        return float(self.hapinitial(a))*float(self.hapinitial(b))

    #Emission state probabilities. gt: genotype. s: tuple of alleles.
    def emission(self, i,s):
        #If one allele is unknown. If known allele is contained in s then 1 is returned.
        if self.gt[i].count('.') == 1:
            for i in self.gt[i]:
                if (i == "|") or (i == "\\"):
                    continue
                if i != '.':
                    if i in s:
                        return 1
                    else:
                        return 0
                
        #If both alleles are unknown 1 is returned
        elif self.gt[i].count('.') == 2:        
            return 1

        #If gt equals s without order 1 is returned
        else:
            if set([self.gt[i][0],self.gt[i][-1]]) == set(s):
                return 1
            else:
                return 0

    #Transition state probabilities. e,d: edge tuples incl data
    def haptrans(self, (e,d)):
        #If parent node of edge e is child node of edge d.
        if e[0] == d[1]:

            #Edge count/parent node count is returned
            return float(e[3]['weight'])/float(sum(self.T[0].node[e[0]]['hap'].itervalues())) 
        else:
            return 0.0

    #Diploid transition probabilities
    def diptrans(self, a,b):
        return float(self.haptrans(a))*float(self.haptrans(b))

    #Function to find edge corresponding to matrix coordinate. g=flattened index. l=level(index of m), e = number of edges
    def findedge(self,g,l,e):    
        #Convert flattened index to coordinate
        t = np.unravel_index(g, (e,e))    
        f = [-1,-1]
        if l == 0:
            x = 1
        else:
            x = [j for j in range(self.T[1][l-1]+1,self.T[1][l]+1)]        
        
        #Set f to equal tuple of ordered edges that the index g corresponds to
        for a, b in enumerate(self.T[0].out_edges(nbunch=x, keys=True, data=True)):
            if t[0] == a:
                f[0] = b
            if t[1] == a:
                f[1] = b

        if -1 in f:
            print 'Error in findedge'
        else:
            return f
    

def treealgorithm(h):      

    #Function which takes node and splits node according to haplotype dictionary attached to the node
    def nodesplit(G,n,l,level):

        #Loop over keys in no,de's haplotype dictionary
        for key, value in G.node[n]['hap'].iteritems():
            #If there exists an edge which corresponds to first character of key.
            for u,v,k,edata in G.out_edges(n,data=True, keys=True):
                if key[0] == edata['a']:
                    #Add value to weight variable on edge
                    G.edge[u][v][k]['weight'] += value
                    #Add haplotype suffix and value to dictionary of connecting node
                    G.node[v]['hap'][key[1:]]=value
                    break

            #If there does not exist an edge.
            else:
                m = level[l+1]+1                              
                #Create edge and label it with first character of key
                G.add_edge(n,m,a=key[0],weight=value)
                #Assign dictionary of suffix of that key to node that has been created.
                G.add_node(m, hap={key[1:]:value})
           
                #Add node to gl
                for i in range(l+1, len(level)):
                    level[i] += 1

    #Functinon tests if two nodes are similar enough to merge. Returns similarity score or false.
    def mergetest(a,b,l):
        #Components of formula are calculated
        ta = sum(G.node[a]['hap'].itervalues())
        tb = sum(G.node[b]['hap'].itervalues())
        thr = math.sqrt(math.pow(ta, -1) + math.pow(tb, -1))
        maxs = 0
        al = [1]*(hlength-l+2)
        bl = [1]*(hlength-l+2)

        #Subgraphs A and B are created in order to calculate similarity score
        A=nx.MultiDiGraph()
        A.clear()
        A.add_node(1,hap=G.node[a]['hap'])
        B=nx.MultiDiGraph()
        B.clear()
        B.add_node(1,hap=G.node[b]['hap'])
    
        q = [[1,1,0]]
    
        while len(q) != 0:

            #Splits nodes subgraphs
            nodesplit(A,q[0][0],q[0][2],al)
            nodesplit(B,q[0][1],q[0][2],bl)

            da = dict((edata['a'],sum(A.node[v]['hap'].itervalues())) for u,v,edata in A.out_edges(q[0][0], data=True))
            db = dict((edata['a'],sum(B.node[v]['hap'].itervalues())) for u,v,edata in B.out_edges(q[0][0], data=True))
    
            sa = set(da.keys())
            sb = set(db.keys())
            
            #Iterates through alelles which are not a member of both node's outgoing edges and tests score against threshold
            for i in list(sa.symmetric_difference(sb)):
            
                if i in da:
                    s = math.fabs(float(da[i])/float(ta))
                    if s > thr:
                         return False
                    elif s > maxs:
                        maxs = s

                if i in db:
                    s = math.fabs(float(db[i])/float(tb))
                    if s > thr:
                        return False
                    elif s > maxs:
                        maxs = s

            #Iterates through alleles which are a member of both node's outgoing edges
            for i in list(sa.intersection(sb)):
                s = math.fabs(float(da[i])/ta - (float(db[i])/tb))
                if s > thr:
                    return False                
                else:
                    if s > maxs:
                        maxs = s              
                    a = [v for u,v,edata in A.out_edges(q[0][0],data=True) if edata['a']==i]                               
                    b = [v for u,v,edata in B.out_edges(q[0][0],data=True) if edata['a']==i]
                
                    #Adds next nodes to list if edges has passed edge test            
                    if len(random.choice(A.node[a[0]]['hap'].keys())) != 0:
                        q.append([a[0],b[0],q[0][2]+1])

            q.pop(0)
          
        return maxs   

    #Function merges 2 nodes. a should always be < b
    def mergenodes(G,a,b):
  
        #Iterate through haplotypes on second node
        for key, value in G.node[b]['hap'].iteritems():
            #If haplotype exists in dictionary of node a. Add weight to dictionary
            if key in G.node[a]['hap']:
                G.node[a]['hap'][key] += value
            #Otherwise add key and value to dictionary of node a
            else:
                G.node[a]['hap'].update({key:value})

        #Move all incoming edges of b to a
        for i in G.in_edges(b, data=True, keys=True):
            G.add_edge(i[0],a,a=i[3]['a'],weight=i[3]['weight'])
            G.remove_edge(i[0],i[1],key=i[2])

        #Remove node b form G
        G.remove_node(b)
        
        #Remove node from gl
        for i in gl:
            if b <= i:
                gl[gl.index(i)] = gl[gl.index(i)] - 1

        nodes = list(set(G.nodes()) - set([i for i in range(b+1)]))        
        mapping=dict(zip(nodes,range(b,len(G.nodes())+1)))
        
        #Relablel nodes so that they are consecutive integers
        G = nx.relabel_nodes(G, mapping, copy=False)


    #Merge function carries out pairwise test between all nodes on given level and merges the lowest scoring nodes.
    #Cycle is repeated until no more merges can be made on the give level.
    def merge(G, l):
        #Set variable nn for number of nodes
        nn = gl[l] - gl[l-1]
    
        #If only one node on level.
        if nn == 1:
            #No need to carry out merge function so returns.
            return
    
        #If there is only 2 nodes in level. 
        elif nn == 2:
            #No need to carry out merge if nodes fail mergetest
            if mergetest(gl[l]-1,gl[l],l) == False:
                return
            #Merge nodes
            else:
                mergenodes(G, gl[l]-1,gl[l])

        #If there are more than 2 nodes in a level, all pairs of nodes are compared and lowest scoring is merged.
        else:
            merge = True
            while merge == True:
                merge = False
                levelmin = 10000000
                minj = 0
                mink = 0
            
                #Calculate merge score for each pair of nodes
                for j in range(gl[l-1]+1, gl[l]+1):
                    for k in range(j, gl[l]+1):                    
                        if j != k:                        
                            testscore = mergetest(j,k,l)                      
                            if testscore == False:
                                continue
                            else:

                                #Details of lowest scoring pair are stored
                                if testscore < levelmin:
                                    levelmin = testscore
                                    minj = j
                                    mink = k                
                                
                #The lowest scoring pair of nodes are merged
                if levelmin != 10000000:
                    
                    mergenodes(G, minj,mink)                
                    #Whenever a merge has taken place, merge is set to True so that the loop is repeated
                    merge = True

    #Create list of first node in each level
    gl=[1]*(hlength+2)
    
    #Create networkx MultiGraph
    G=nx.MultiDiGraph()

    #Add start node 1 and split the first node.
    G.add_node(1,hap=h)
    nodesplit(G,1,0,gl)
    merge(G,1)

    for i in range(1,hlength):

        for j in range(gl[i-1]+1,gl[i]+1):
            nodesplit(G,j,i,gl)
       
        merge(G,i+1)

    #Set endnode as first node in last level
    endnode = gl[-1]+1

    #Iterate through all nodes on penultimate level
    for i in range(gl[-3]+1,gl[-2]+1):
        for key, value in G.node[i]['hap'].iteritems():
            G.add_edge(i,endnode,a=key,weight=value)

    #Set glevel so it is correct    
    gl[-1] = endnode    
    return (G, gl)


def forwardbackward(T, gt):
    fb = HMM(T, gt)

    #Create n the list of the number of edges in each level
    n = [T[0].out_degree(1)]

    #m is the list of matrices of forward probabilities at each level
    m = [np.zeros(shape=(n[0],n[0]))]

    #Append matrices to list m
    for i in range(hlength):    
        #Sum n for each level
        n.append(T[0].out_degree(T[1][i]+1))
        for j in range(T[1][i]+2, T[1][i+1]+1):
            n[i+1] += T[0].out_degree(j)

        #Append zero-filled matrix to list m
        m.append(np.zeros(shape=(n[i+1],n[i+1])))
 
    #Initiation. Iterate through pairs of outgoing edges from node 1.
    for a, b  in itertools.product([(i, j) for i, j in enumerate(T[0].out_edges(1, keys=True, data=True))], repeat=2):

        #If emmsion probability does not equal 0
        if fb.emission(0,(a[1][3]['a'],b[1][3]['a'])) != 0:              
            if a[1] == b[1]:
                t = fb.hapinitial(a[1][3]['a'])
                var = t*t            
            else:
                var = fb.dipinitial(a[1][3]['a'], b[1][3]['a'])

            #Matrix element is set to calculated diploid initial probability
            m[0][a[0]][b[0]] = var


    #Induction. Iterate through each level
    for i in range(1,hlength+1):
        #Iterate through ordered pairs of outgoing edges in each level
        for a, b in itertools.product([(c, d) for c, d in enumerate(T[0].out_edges(nbunch=[j for j in range(T[1][i-1]+1,T[1][i]+1)], keys=True, data=True))], repeat=2):


            #If emission probability does not equal 0
            if fb.emission(i,(a[1][3]['a'],b[1][3]['a'])) != 0.0:
                if i == 1:
                    nodes = 1
                else:                
                    nodes = [j for j in range(T[1][i-2]+1,T[1][i-1]+1)]

                var = sum([(m[i-1][c[0]][d[0]]*fb.diptrans([a[1],c[1]], [b[1],d[1]])) for c, d  in itertools.product([(e, f) for e, f in enumerate(T[0].out_edges(nbunch=nodes, keys=True, data=True))], repeat=2)])
               
                #Matrix element is set to var calculation formula
                m[i][a[0]][b[0]] = var


    #Backwards sampling
    #Initial probability
    inprob = []

    #Create probability distribution for sampling first allele
    for i in m[3].flat:
        inprob.append(i/m[hlength].sum())

    #List of  flattened indices chosen from sampling 
    s = [0]*(hlength+1)

    #Corresponding probabilities of the chosen positions 
    p = [0]*(hlength+1)


    #Randomly choose first element in m[3] according to initial probabilities
    s[hlength] = np.random.choice(m[hlength].size, p=inprob)

    #Set probability of the chosen position
    p[hlength] = inprob[s[hlength]]

    #Iterate through levels in reverse order
    for i in range(hlength, 0, -1):

        prob = []
        #Set e to equal the edge description of edges sampled

        e = fb.findedge(s[i], i, n[i])

        #For each position in forward probability matrix of level below, calculate sampling probabilities
        for b, j in enumerate(m[i-1].flat):
            d = fb.findedge(b, i-1, n[i-1])
            prob.append((fb.emission(i,(e[0][3]['a'],e[1][3]['a']))*fb.diptrans((e[0],d[0]),(e[1],d[1]))*j)/m[i].flat[s[i]])

        #Choose next sampled edge based on calculated probabilities
        s[i-1] = np.random.choice(m[i-1].size, p=prob)

        #Record probability of chosen edges.
        p[i-1] = prob[s[i-1]]

    #Print descriptive list of chosen edges
    #for i in range(hlength+1):
    #print findedge(s[i], i)
    #print 'The pair of paths has sampling probability in either order of'

    #Calculate probability of sampled path
    #print np.product(p)*2

    sample = ['','']

    for i in range(hlength+1):
        e = fb.findedge(s[i],i, n[i])
        sample[0] += e[0][3]['a']
        sample[1] += e[1][3]['a']

    return sample

def viterbi(T, gt):
    vi = HMM(T, gt)
    #Create n the list of the number of edges in each level
    n = [T[0].out_degree(1)]
    
    #v is the list of matrices of viterbi probabilities at   m each level
    v = [np.zeros(shape=(n[0],n[0]))]


    #arglist is the list of matrices of path to find maximum viterbi probabilities
    arglist = [np.zeros((n[0],n[0]))]
    arglist[0].fill(-1)

    #Append matrices to list arglist
    for i in range(hlength):
        #Sum n for each level
        n.append(T[0].out_degree(T[1][i]+1))
        for j in range(T[1][i]+2, T[1][i+1]+1):
            n[i+1] += T[0].out_degree(j)
            
        v.append(np.zeros(shape=(n[i+1],n[i+1])))
        arglist.append(np.zeros((n[i+1],n[i+1])))
        arglist[i+1].fill(-1)

    #Initiation. Iterate through pairs of outgoing edges from node 1. SAME AS FORWARD ALGORITHM.
    for a, b in itertools.product([(i, j) for i, j in enumerate(T[0].out_edges(1, keys=True, data=True))], repeat=2):
    
        #If emmsion probability does not equal 0
        if vi.emission(0,(a[1][3]['a'],b[1][3]['a'])) != 0:               
            if a[1] == b[1]:
                t = vi.hapinitial(a[1][3]['a'])
                var = t*t            
            else:
                var = vi.dipinitial(a[1][3]['a'], b[1][3]['a'])

            #Matrix element is set to calculated diploid initial probability
            v[0][a[0]][b[0]] = var


    #Induction. Iterate through each level
    for i in range(1,hlength+1):

        #Iterate through ordered pairs of outgoing edges in each level
        for a, b in itertools.product([(c, d) for c, d in enumerate(T[0].out_edges(nbunch=[j for j in range(T[1][i-1]+1,T[1][i]+1)], keys=True, data=True))], repeat=2):       

            #If emission probability does not equal 0
            if vi.emission(i,(a[1][3]['a'],b[1][3]['a'])) != 0:
            
                if i == 1:
                    nodes = 1
                else:
                    nodes = [j for j in range(T[1][i-2]+1,T[1][i-1]+1)]

                max = 0
                args = -1

                #Calculate maximum value and record which edges correspond.
                for c, d  in itertools.product([(e, f) for e, f in enumerate(T[0].out_edges(nbunch=nodes, keys=True, data=True))], repeat=2):                    
                    value = (v[i-1][c[0]][d[0]]*vi.diptrans([a[1],c[1]], [b[1],d[1]]))
                   
                    if max < value:
                        max = value                    
                        args = np.ravel_multi_index((c[0],d[0]), (n[i-1],n[i-1]), mode='raise')                   
                  
                #Matrix element is set to maximum and element this comes from is recorded.
                v[i][a[0]][b[0]] = max
                
                arglist[i][a[0]][b[0]] = args


    #Backtracking process
    value =  v[hlength].argmax()

    phased = []

    #Iterate through levels bacwards and create list of phased alleles.
    for i in range(hlength,-1, -1):
        edge = vi.findedge(value, i, n[i])
        phased.insert(0,(edge[0][3]['a'],(edge[1][3]['a'])))
        value = int(arglist[i].flat[value])


    return phased
    
    

def reverseorder(haplotypes):
    reversehaplotype = {}
    for key, value in haplotypes.iteritems():
        reversehaplotype.update({key[::-1] : value})
    return reversehaplotype
    
def treesequence(haplotypes,r):
    #Input into tree algorithm
    T = treealgorithm(haplotypes)
    
    haplotypes={}
    for i in GT:        
        if r == 1:
            i = i[::-1]
                    
        samples = forwardbackward(T,i)

        for j in samples:

            if j in haplotypes:
                haplotypes[j] += 1
            else:
                haplotypes[j] = 1
    return haplotypes


#Extract genotypes from data
with open('/Users/mp18/Documents/bloch/Data/Table_1', 'rb') as f:
    reader = csv.reader(f, delimiter='\t',skipinitialspace = True)
    #Create list of tuples which contain genotypes of each individual
    GT = zip(*reader)

#Remove first tuple of position names
del GT[0]

#Create input for tree algorithm
haplotypes = {}

#Randomly assign phase for each sample
for i in GT:
    
    a=''
    b=''
    for j in i:
        x=j[0]
        y=j[-1]
        rand = random.randrange(0,2)
        if rand == 0:
            a+=x
            b+=y
        else:
            a+=y
            b+=x
    
    if a in haplotypes:
        haplotypes[a] += 1
    else:
        haplotypes[a] = 1

    if b in haplotypes:
        haplotypes[b] += 1
    else:
        haplotypes[b] = 1


#Haplotype length
hlength = len(a)-1
print "initial haplotypes"

haplotypes = {'0110': 53, '0111': 25, '0000': 27, '0001': 70, '0011': 96, '0010': 28, '0101': 6, '0100': 6, '1111': 9, '1110': 10, '1100': 4, '1101': 3, '1010': 15, '1011': 127, '1001': 102, '1000': 19}
print haplotypes
haplotypes = treesequence(haplotypes, 0)

print "first iteration"
print haplotypes

iterations = 1

r = 1
#n must be odd number so that final iteration haplotypes are correct way round.
m=1

while iterations < m:
    print iterations
    print "haplotypes in"
    print haplotypes
    haplotypes = reverseorder(haplotypes)
    print "haplotypes after reverse"
    print haplotypes
    print "r"
    print r
    haplotypes = treesequence(haplotypes,r)
    print "haplotypes out"
    print haplotypes
    iterations += 1
    r += 1
    r = r%2
   
print "last haplotypes"
print haplotypes

#Input into tree algorithm
T = treealgorithm(haplotypes)
    
phased=[]
for i in GT:
    print i
    result = viterbi(T, i)
    print result
    for j in result:
        pass
