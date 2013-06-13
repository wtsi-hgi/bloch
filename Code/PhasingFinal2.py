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
import sys


class HMM:
    def __init__(self, T, gt):
        self.T = T
        self.gt = gt
        
    #Haploid initial state probabilities. allele: allele number.
    def hapinitial(self, allele):
        #Edge counts are used. Count of edge/Total count for all edges.
        for i in self.T[0].out_edges(1, data=True):        
            if i[2]['a'] == allele:                
                return float(i[2]['weight'])/float(self.T[0].node[1]['weight'])

    #Diploid initial state probabilities. a,b: allele number
    def dipinitial(self, a,b):    
        return float(self.hapinitial(a))*float(self.hapinitial(b))

    #Emission state probabilities. gt: genotype. i: level. s: tuple of alleles.
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
            return float(e[3]['weight'])/float(self.T[0].node[e[0]]['weight']) 
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
        #if l == 0:
        #    x = 1
        #else:
        #    x = T[1][l]        
 
        #Set f to equal tuple of ordered edges that the index g corresponds to
        for a, b in enumerate(self.T[0].out_edges(nbunch=self.T[1][l], keys=True, data=True)):
          
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
    def nodesplit(G,n,nodes):

        l = G.node[n]['level']
        
        #Loop over keys in node's haplotype dictionary
        for key, value in G.node[n]['hap'].iteritems():
            #If there exists an edge which corresponds to first character of key.
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
                #Assign label of new node
                if len(nodes[l+1]) != 0:
                    m = max(nodes[l+1]) + 1
                else:
                    m = max(nodes[l]) + 1
                #Create edge and label it with first character of key
                G.add_edge(n,m,a=key[0],weight=value)
                #Assign dictionary of suffix and level marker of that key to node that has been created.
                G.add_node(m, hap={key[1:]:value},level=l+1)
                #Add node to list of nodes
                nodes[l+1].append(m)
        
        G.add_node(n, weight=sum(G.node[n]['hap'].values()))
        #Remove haplotype information of node that has been split
        del G.node[n]['hap']


    #Functinon tests if two nodes are similar enough to merge. Returns similarity score or false.
    def mergetest(a,b):
        
        if G.node[a]['level'] != G.node[b]['level']:
            raise ValueError('Nodes to tested are not on the same level')
        
        #print a, b
        #print G.node[a], G.node[b]
        #Components of formula are calculated
        ta = sum(G.node[a]['hap'].itervalues())
        tb = sum(G.node[b]['hap'].itervalues())
        thr = math.sqrt(math.pow(ta, -1) + math.pow(tb, -1))
        maxs = 0
        anodes = []
        bnodes = []

        for i in range(hlength-G.node[a]['level']+2):
            anodes.append([])
            bnodes.append([])

        #Subgraphs A and B are created in order to calculate similarity score
        A=nx.MultiDiGraph()
        A.clear()
        A.add_node(1,hap=G.node[a]['hap'],level=0)
        anodes[0].append(1)
        B=nx.MultiDiGraph()
        B.clear()
        B.add_node(1,hap=G.node[b]['hap'],level=0)
        bnodes[0].append(1)
    
        q = [[1,1]]
    
        while len(q) != 0:
            #print q
            #Splits nodes subgraphs
            nodesplit(A,q[0][0],anodes)
            nodesplit(B,q[0][1],bnodes)

            #print "Graph A"
            #for i in A.nodes(data=True):
            #    print i
            #for i in A.edges(data=True, keys=True):
            #    print i
            #print "Graph B"
            #for i in B.nodes(data=True):
            #    print i
            #for i in B.edges(data=True, keys=True):
            #    print i                                
            
            da = dict((edata['a'],sum(A.node[v]['hap'].itervalues())) for u,v,edata in A.out_edges(q[0][0], data=True))
            db = dict((edata['a'],sum(B.node[v]['hap'].itervalues())) for u,v,edata in B.out_edges(q[0][1], data=True))          
            
            sa = set(da.keys())           
            sb = set(db.keys())            
            
            #Iterates through alelles which are not a member of both node's outgoing edges and tests score against threshold
            for i in sa.symmetric_difference(sb):                
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
            for i in sa.intersection(sb):
                s = math.fabs(float(da[i])/ta - (float(db[i])/tb))
                
                if s > thr:
                    return False                
                else:
                    if s > maxs:
                        maxs = s
                        
                    a = [v for u,v,edata in A.out_edges(q[0][0],data=True) if edata['a']==i]                               
                    b = [v for u,v,edata in B.out_edges(q[0][1],data=True) if edata['a']==i]
                    
                    #Adds next nodes to list if edges has passed edge test            
                    if len(random.choice(A.node[a[0]]['hap'].keys())) != 0:                        
                        q.append([a[0],b[0]])

            q.pop(0)
          
        return maxs


    #Function merges 2 nodes. a should always be < b
    def mergenodes(a,b):
        if G.node[a]['level'] != G.node[b]['level']:
            raise ValueError('Nodes to merged are not on the same level')        
  
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

            #G.add_node(a, weight=sum(G.node[a]['hap'].values()))
        #Remove node b from G
        G.remove_node(b)
        
        #Remove node b from list of nodes
        gnodes[G.node[a]['level']].remove(b)
        

    #Merge function carries out pairwise test between all nodes on given level and merges the lowest scoring nodes.
    #Cycle is repeated until no more merges can be made on the give level.
    def merge(l):
   
        #If only one node on level.
        if len(gnodes[l]) == 1:
            #No need to carry out merge function so returns.
            return
    
        #If there is only 2 nodes in level. 
        elif len(gnodes[l]) == 2:
            #No need to carry out merge if nodes fail mergetest
            if mergetest(gnodes[l][0],gnodes[l][1]) == False:
                return
            #Merge nodes
            else:
                mergenodes(gnodes[l][0],gnodes[l][1])

        #If there are more than 2 nodes in a level, all pairs of nodes are compared and lowest scoring is merged.
        else:
            merge = True
            while merge == True:
                merge = False
                levelmin = 10000000
                minj = 0
                mink = 0

                for j, k in itertools.combinations(gnodes[l], 2):                  
                    
                    testscore = mergetest(j,k)                      
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
                    
                    mergenodes(minj,mink)

                    #Whenever a merge has taken place, merge is set to True so that the loop is repeated
                    merge = True
    
        
    gnodes = []
    #Create list of nodes on each level.
    for i in range(hlength+2):
        gnodes.append([])    

    #Create networkx MultiGraph
    G=nx.MultiDiGraph()
    
    #Add start node 1 and split the first node.
    G.add_node(1,hap=h,level=0)
    
    gnodes[0].append(1)    
    
    nodesplit(G,1,gnodes)

    merge(1)

    for i in range(1,hlength):        
        for j in gnodes[i]:            
            nodesplit(G,j,gnodes)        
        merge(i+1)

    #Set endnode as first node in last level
    endnode = max(gnodes[-2]) +1

    #Iterate through all nodes on penultimate level
    for i in gnodes[-2]:
        for key, value in G.node[i]['hap'].iteritems():
            G.add_edge(i,endnode,a=key,weight=value)
        G.add_node(i, weight=sum(G.node[i]['hap'].values()))
        del G.node[i]['hap']

    G.node[endnode]['level'] = hlength+1

    #Set gnodes so it is correct    
    gnodes[-1].append(endnode)
    for i in G.nodes(data=True):
        print i

    for i in G.edges(data=True, keys=True):
        print i

    return (G, gnodes)

def forwardbackward(T, gt):

    def random_weighted_choice(weights):
        rnd = random.random() * sum(weights)
        for i, w in enumerate(weights):
            rnd -= w
            if rnd < 0:
                return i
        
    fb = HMM(T, gt)

    #Create n the list of the number of edges in each level
    n = [T[0].out_degree(1)]


    #m is the list of matrices of forward probabilities at each level
    m = [np.zeros(shape=(n[0],n[0]))]

    #Append matrices to list m
    for i in range(hlength):    
        #Sum n for each level
        n.append(0)
        for j in T[1][i+1]:
            n[i+1] += T[0].out_degree(j)
        
        #Append zero-filled matrix to list m
        m.append(np.zeros(shape=(n[i+1],n[i+1])))

    
    #Initiation. Iterate through pairs of outgoing edges from node 1.
    for a, b in itertools.product([(i, j) for i, j in enumerate(T[0].out_edges(1, keys=True, data=True))], repeat=2):
        
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
        for a, b in itertools.product([(c, d) for c, d in enumerate(T[0].out_edges(nbunch=T[1][i], keys=True, data=True))], repeat=2):
            
            #If emission probability does not equal 0
            if fb.emission(i,(a[1][3]['a'],b[1][3]['a'])) != 0.0:

                var = sum([(m[i-1][c[0]][d[0]]*fb.diptrans([a[1],c[1]], [b[1],d[1]])) for c, d  in itertools.product([(e, f) for e, f in enumerate(T[0].out_edges(nbunch=T[1][i-1], keys=True, data=True))], repeat=2)])
               
                #Matrix element is set to var calculation formula
                m[i][a[0]][b[0]] = var


    #Backwards sampling


    #List of  flattened indices chosen from sampling 
    s = [0]*(hlength+1)

    #Corresponding probabilities of the chosen positions 
    p = [0]*(hlength+1)
    #Randomly choose first element in m[3] according to initial probabilities
    s[hlength] = random_weighted_choice(m[hlength].flatten())
    

    #Set probability of the chosen position
    p[hlength] = m[hlength].flatten()[s[hlength]]

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
        s[i-1] = random_weighted_choice(prob)

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
    n = [T[0].node[1]['weight']]
    
    #v is the list of matrices of viterbi probabilities at   m each level
    v = [np.zeros(shape=(n[0],n[0]))]


    #arglist is the list of matrices of path to find maximum viterbi probabilities
    arglist = [np.zeros((n[0],n[0]))]
    arglist[0].fill(-1)

    #Append matrices to list arglist
    for i in range(hlength):
        #Sum n for each level
        n.append(0)
        for j in T[1][i+1]:
            n[i+1] += T[0].node[j]['weight']
            
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
        for a, b in itertools.product([(c, d) for c, d in enumerate(T[0].out_edges(nbunch=T[1][i], keys=True, data=True))], repeat=2):       

            #If emission probability does not equal 0
            if vi.emission(i,(a[1][3]['a'],b[1][3]['a'])) != 0:

                max = 0
                args = -1

                #Calculate maximum value and record which edges correspond.
                for c, d  in itertools.product([(e, f) for e, f in enumerate(T[0].out_edges(nbunch=T[1][i-1], keys=True, data=True))], repeat=2):                    
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

f.close()
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
iterations = 1
r=1
#m has to be an odd number
m=1

haplotypes = treesequence(haplotypes, 0)
sys.stdout.write(str(m)+" iterations altogether\n1\n")


while iterations < m:
    sys.stdout.write(str(iterations+1)+"\n")    
    
    haplotypes = reverseorder(haplotypes)
    haplotypes = treesequence(haplotypes,r)
    iterations += 1
    r += 1
    r = r%2

    #print "last haplotypes"
    #print haplotypes

#Input into tree algorithm
T = treealgorithm(haplotypes)
    
phased=[]
correct = 0
incorrect = 0

for i in GT:
    result = viterbi(T, i)
    
    for j in range(len(result)):
        
        if result[j][0] == i[j][0]:
            correct += 1
        else:
            incorrect += 1
            
        if result[j][1] == i[j][-1]:
            correct += 1
        else:
            incorrect += 1

g = open("output.txt", "w")
g.write("Correct: "+ str(correct)+"\n")
g.write("Incorrect: "+str(incorrect))

g.close()
