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
import pickle
import cPickle
import string
import functools
import argparse
import time

#Funciton which performs random weighted choice
def random_weighted_choice(weights):
    rnd = random.random() * sum(weights)
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return i
        
#Tree class. T = Tree. tnodes = list of nodes. d = direction. 0 = forward. 1 = backwards. c = node count.
class Tree:
    def __init__(self, direction):
        #Initialize variables
        self.T = nx.MultiDiGraph()
        self.tnodes = []
        self.d = direction
        self.c = 1        
        
        #Add first node
        self.T.add_node(1,level=0,weight=0)
        
        self.tnodes = [[] for x in xrange(hlength+2)]
#        for i in range(hlength+2):
#            self.tnodes.append([])
        self.tnodes[0].append(1)
        
    #Method which adds a genotype to the graph. Input of first iteration has already been randomly phased and imputed.
    #gt is a list of lists. Each marker corresponds to a list with 2 elements in.
    def add_genotype(self, gt):
        
        #Add single marker position to the graph. cnode = current node. a = allele. l = level
        def add_marker(cnode,a,l):
            #If current node has no outgoing edges
            if self.T.out_degree(cnode) == 0:
                #Increase node name counter by one
                self.c += 1
                
                #Add edge between current node and new node.
                self.T.add_edge(cnode,self.c,a=a,weight=1)
                self.T.add_node(self.c, level=l+1,weight=1)
                self.tnodes[l+1].append(self.c)
                
                #Return current node position
                return self.c
            
            else:
                #Search through all outgoing nodes to find corresponding allele
                for u,v,k,edata in self.T.out_edges(cnode,data=True,keys=True):     
                    #If corresponding allele found. Add weight to edge and node count.           
                    if edata['a'] == a:
                        self.T.edge[u][v][k]['weight'] += 1
                        self.T.node[v]['weight'] +=1
                        
                        #Return current node position
                        return v
                
                #If no match is found add new edge and node
                self.c +=1
                self.T.add_edge(cnode,self.c,a=a,weight=1)
                self.T.add_node(self.c,level=l+1,weight=1)
                self.tnodes[l+1].append(self.c)
                
                #Return current node position
                return self.c               
            
        #Weight added to first node
        self.T.node[1]['weight'] += 2   
             
        #Current nodes set for haplotypes a and b
        anode = 1
        bnode = 1
        
        #Iterate through each marker in the position and add to the graph
        for n, i in enumerate(gt):        
            anode = add_marker(anode, i[0], n)           
            bnode = add_marker(bnode, i[1], n)
        
    #Merge nodes on each level
    def merge(self):

        def mergetest(a,b):
            #Components of formula are calculated
            ta = self.T.node[a]['weight']
            tb = self.T.node[b]['weight']
            thr = math.sqrt(math.pow(ta, -1) + math.pow(tb, -1))
            maxs = 0        
            
            q=[[a,b]]       
    
            while len(q) != 0:
    
                da = dict((edata['a'],self.T.edge[u][v][k]['weight']) for u,v,k,edata in self.T.out_edges(q[0][0], data=True, keys= True))
                db = dict((edata['a'],self.T.edge[u][v][k]['weight']) for u,v,k,edata in self.T.out_edges(q[0][1], data=True, keys = True))    
    
                sa = set(da.keys())           
                sb = set(db.keys())         
    
                #Iterate through elements that are in sa that are not in sb
                for i in sa.difference(sb):                     
                    s = math.fabs(float(da[i])/float(ta))                
                    if s > thr:
                        return 0.0
                    elif s > maxs:
                        maxs = s
    
                #Iterate through elements that are in sb but not sa
                for i in sb.difference(sa):        
                    s = math.fabs(float(db[i])/float(tb))                
                    if s > thr:
                        return 0.0
                    elif s > maxs:
                        maxs = s
    
                #Iterates through alleles which are a member of both node's outgoing edges
                for i in sa.intersection(sb):
                    s = math.fabs(float(da[i])/ta - (float(db[i])/tb))
                    
                    if s > thr:
                        return 0.0             
                    else:
                        if s > maxs:
                            maxs = s                        
                        a = [v for u,v,edata in self.T.out_edges(q[0][0],data=True) if edata['a']==i]                   
                        b = [v for u,v,edata in self.T.out_edges(q[0][1],data=True) if edata['a']==i]
                        
                        #Adds next nodes to list if edges has passed edge test            
                        if self.T.node[a[0]]['level'] != (hlength+1):                        
                            q.append([a[0],b[0]])
    
                q.pop(0)
            return maxs
    
        def mergenodes(a,b):
            #Move all incoming edges of b to a
            for i,j,k,edata in self.T.in_edges(b, data=True, keys=True):
                self.T.add_edge(i,a,a=edata['a'],weight=edata['weight'])
                self.T.node[a]['weight'] += self.T.edge[i][j][k]['weight']
                self.T.remove_edge(i,j,key=k)
    
            q=[[a,b]]    
            while len(q) != 0:    
                #Pairwise iterate through outgoing edges of nodes at front of queue
                for ba, bb, bk, bedata in self.T.out_edges(q[0][1], data=True, keys=True):    
                    found = 0
                    for aa, ab, ak, aedata in self.T.out_edges(q[0][0], data=True, keys=True):
                        #If edge exists on both nodes    
                        if bedata['a'] == aedata['a']:
                            found += 1
                            self.T.edge[aa][ab][ak]['weight'] += bedata['weight']
                            self.T.node[ab]['weight'] += bedata['weight']
                            if self.T.node[bb]['level'] != (hlength+2):  
                                q.append([ab,bb])                                
                            break               
                    #If edge does not exist
                    if found == 0:
                        self.T.add_edge(q[0][0],bb,a=bedata['a'],weight=bedata['weight'])
                        
                    #Remove edge which has been merged    
                    self.T.remove_edge(ba,bb,key=bk)
                    if self.T.node[bb]['level'] == (hlength+2):
                        self.tnodes[self.T.node[bb]['level']].remove(bb)
                        self.T.remove_node(bb)
    
                self.tnodes[self.T.node[q[0][1]]['level']].remove(q[0][1])
                self.T.remove_node(q[0][1])
                q.pop(0)
                
            
                
        
            
        #Iterate through each level               
        for i in xrange(hlength):
            #No merge on level if only one node
            if len(self.tnodes[i]) == 1:
                pass
            #One mergetest required if only two nodes on level
            elif len(self.tnodes[i]) == 2:
                if mergetest(self.tnodes[i][0],self.tnodes[i][1]) != 0.0:
                    mergenodes(self.tnodes[i][0],self.tnodes[i][1])
            #Node similarity scores are stored in a matrix. Lowest scoring pair is merged and
            #scores recalculated.
            else:
                #Create copy of list of node names on the level
                nlist = list(self.tnodes[i])                
                #n is the number of nodes on that level
                n = len(nlist)    
                #k list of 0s to mark which nodes have been deleted
                k = np.zeros(n) 
                #Create empty numpy matrix of similarity scores
                simmatrix = np.zeros((n,n))
    
                #Calculate values in half of matrix
                for a, b in itertools.combinations(enumerate(self.tnodes[i]), 2):
                    simmatrix[a[0],b[0]] = mergetest(a[1],b[1])                
                 
                #While some values in matrix are not equal to zero
                while np.count_nonzero(simmatrix) != 0:
                    #Mask 0s to find minimum value which is not 0
                    ma = np.ma.masked_equal(simmatrix, 0.0, copy=False)
                    #Return 2D coordinates for position of minimum score
                    mnodes = np.unravel_index(ma.argmin(), (n,n))        
                    #Merge nodes with lowest score
                    mergenodes(nlist[mnodes[0]],nlist[mnodes[1]])
                    #deleted nodes need to be marked
                    k[mnodes[1]] = 1
    
                    #Set all matrix positions of deleted node to 0
                    for j in xrange(len(k)):
                        if j <= mnodes[1]:
                            simmatrix[j][mnodes[1]] = 0
                        elif j > mnodes[1]:
                            simmatrix[mnodes[1]][j] = 0
                                          
                    #Recalculate similarity score for changed node
                    for j in xrange(len(k)):
                        if k[j] == 1:
                            continue
                        else:
                            if j == mnodes[0]:
                                continue
                            elif j < mnodes[0]:
                                simmatrix[j][mnodes[0]] = mergetest(nlist[j],nlist[mnodes[0]])
                            elif j > mnodes[0]:
                                simmatrix[mnodes[0]][j] = mergetest(nlist[mnodes[0]],nlist[j])
                        
        
        
        #Merge all node on the final level   
        self.c += 1
        self.T.add_node(self.c, weight = 0, level=hlength+1)
    
        for i in self.tnodes[-1]:
            
            for a,b,k,edata in self.T.in_edges(i, data=True, keys=True):
                
                self.T.add_edge(a,self.c,a=edata['a'],weight=edata['weight'])
                
                self.T.node[self.c]['weight'] += edata['weight']
                
                self.T.remove_edge(a,b,key=k)
                
            
            self.T.remove_node(i)
        
        self.tnodes[-1] = self.c                         
                
        return self
        
        
class HMM:
    def __init__(self,Tree,gt):
        self.T = Tree.T
        self.gt = gt
    

    #Diploid initial state probabilities. a,b: allele number
    def dipinitial(self, a,b):
        for i in self.T.out_edges(1,data=True):
            if i[2]['a'] == a:
                valuea = float(i[2]['weight'])/float(self.T.node[1]['weight'])
            if i[2]['a'] == b:
                valueb = float(i[2]['weight'])/float(self.T.node[1]['weight'])
                
        return valuea*valueb        
 

    #Emission state probabilities. i: level. s: tuple of alleles.
    def emission(self, i,s):
        genotype = [self.gt[i][0], self.gt[i][-1]]
        genotype.sort()        
        
        if genotype[0] == '.':
            if genotype[1] == '.':
                return 1.0
            else:
                if genotype[1] in s:
                    return 1.0      
        else:
            if set(genotype) == set(s):
                return 1.0
            else:
                return 0.0
  

    #Diploid transition probabilities
    def diptrans(self, a, b):        
        (p,c,w) = a
        (q,d,x) = b
        if p == c and q == d:
            return (float(w)*float(x))/(float(self.T.node[p]['weight'])*float(self.T.node[q]['weight']))
        return 0.0
   
    
        
def forwardbackward(Tree,gt):
    fb = HMM(Tree,gt)     
    
    #Create empty list containing lists of forward probabilities and edge pairs at each level
    edges=[[]]
    forward = [[]]    
        
    #Initiation. Iterate through pairs of outgoing edges from node 1.
    for a, b in itertools.product([j for j in Tree.T.out_edges(1, keys=True, data=True)], repeat=2):
        var = 0.0
        #If emmsion probability does not equal 0
        if fb.emission(0,(a[3]['a'],b[3]['a'])) != 0.0:              

            var = fb.dipinitial(a[3]['a'], b[3]['a'])
            
        if var != 0.0:
            edges[0].append((a,b))
            forward[0].append(var)    

    #Induction.  Iterate through each level.
    for i in xrange(1,hlength+1):

        edges.append([])
        forward.append([])
      
        for a, b in itertools.product([j for j in Tree.T.out_edges(nbunch=Tree.tnodes[i], keys=True, data=True)], repeat=2):
            var = 0.0
            #If emission probability does not equal 0
            if fb.emission(i,(a[3]['a'],b[3]['a'])) != 0.0:
                
                for j, k in enumerate(edges[i-1]):

                    var += (forward[i-1][j]*(fb.diptrans((a[0],k[0][1],a[3]['weight']),(b[0],k[1][1],b[3]['weight']))))
                    
            if var != 0.0:
                edges[i].append((a,b))
                forward[i].append(var)  
   
           
    #s stores the list of chosen edge pairs from front to back             
    s = []
    #Backward stores the list of probabilities of chosen pairs from front to back as chosen
    backward = []
    
    #Choose first s from last matrix of forward probabilities 
    index = random_weighted_choice(forward[hlength])
    
    back_append = backward.append
    s_append = s.append
     
    back_append(forward[hlength][index])
    s_append(edges[hlength][index])    

    for i in xrange(hlength-1,-1,-1):
        #Temporary list containing probabilities that choice is made from
        prob = []
        prob_append = prob.append
        for j, k in enumerate(edges[i]):           
            prob_append(fb.diptrans((s[hlength-i-1][0][0],k[0][1],s[hlength-i-1][0][3]['weight']),(s[hlength-i-1][1][0],k[1][1],s[hlength-i-1][1][3]['weight']))*forward[i][j]/(backward[hlength-i-1]))
        index = random_weighted_choice(prob)
        back_append(forward[i][index])
        s_append(edges[i][index])
    
    #Create output of path that has been sampled
#    alleles = []
#    for i in s:
#        alleles.append([i[0][3]['a'],i[1][3]['a']])

    alleles = [(i[0][3]['a'],i[1][3]['a']) for i in s]
        
    return alleles
    
def viterbi(Tree, gt):
    vi = HMM(Tree,gt)     
    
    #Create empty list containing lists of forward probabilities and edge pairs at each level
    edges=[[]]
    viterbi = [[]]    
    maximum = [[]]
        
    #Initiation. Iterate through pairs of outgoing edges from node 1.
    for a, b in itertools.product([j for j in Tree.T.out_edges(1, keys=True, data=True)], repeat=2):
        var = 0.0
        #If emmsion probability does not equal 0
        if vi.emission(0,(a[3]['a'],b[3]['a'])) != 0.0:              
            var = vi.dipinitial(a[3]['a'], b[3]['a'])
            
        if var != 0.0:
            edges[0].append((a,b))
            viterbi[0].append(var)    

    #Induction.  Iterate through each level.
    for i in xrange(1,hlength+1):

        edges.append([])
        viterbi.append([])
        maximum.append([])
      
        for a, b in itertools.product([j for j in Tree.T.out_edges(nbunch=Tree.tnodes[i], keys=True, data=True)], repeat=2):
            m = 0.0
            medge = 0
            #If emission probability does not equal 0
            if vi.emission(i,(a[3]['a'],b[3]['a'])) != 0.0:
                
                for j, k in enumerate(edges[i-1]):

                    var = (viterbi[i-1][j]*(vi.diptrans((a[0],k[0][1],a[3]['weight']),(b[0],k[1][1],b[3]['weight']))))
                    if m < var:
                        m = var
                        medge = k
                        
            if m != 0.0:
                edges[i].append((a,b))
                viterbi[i].append(m)
                maximum[i].append(medge)          


    value = viterbi[hlength].index(max(viterbi[hlength]))

    final = []
    
    for i in xrange(hlength,-1,-1):    
        e = edges[i][value]
        final.insert(0,[e[0][3]['a'],e[1][3]['a']])
        
        if i == 0:
            break
        value = edges[i-1].index(maximum[i][value])      
    
    return final
        
    
    
#Beginning of code.  Print start time.

def testroute(G, i):
    first_node = 1
    
#    ht1 = []
#    ht2 = []
#    for pos in i:
#        ht1.append(pos[0])
#        ht2.append(pos[-1]) 
    
    ht1 = [pos[0] for pos in i]
    ht2 = [pos[-1] for pos in i]


    def testroute_1(node, ht):
        if ht == []:
            return True
        else:
            next_nodes = [i[1] for i in G.T.out_edges(node, data=True, keys=True) if i[3]['a']== ht[0]]
            if next_nodes == []:
                return False
            for n in next_nodes:
                if testroute_1(n, ht[1:]):
                    return True
            return False
        
    return testroute_1(first_node, ht1) and testroute_1(first_node, ht2)


sys.stdout.write("Start\t"+str(time.clock())+"\n")

#Parse command line arguements
parser = argparse.ArgumentParser()
parser.add_argument('-t','--tree',dest='tree_input',help='input genotype file. phased or unphased.')
parser.add_argument('-o','--output',dest='output',help='output filename')
args = parser.parse_args()   

filename = args.tree_input

#args.tree_input
if filename != None:
    #Create empty list of allele frequencies
    allelefreq=[]
    alleles=[]

    #Extract genotypes from data
    f = open(filename, 'rb')
    #Create marker name list for use at the end
    imputed = []
    for line in f:
        dic = {}   
        l,r= line.split('\t',1)
        imputed.append(l+"\t")
        for i in r:
            if i in ('|','.','/','\t','\n'):
                continue
            else:
                if i in dic:
                    dic[i] += 1
                else:
                    dic[i] = 1
        alleles.append(dic.keys())
        allelefreq.append(dic.values())
    
    #Move current position to beginning of the file
    f.seek(0,0)
    reader = csv.reader(f, delimiter='\t',skipinitialspace = True)        
    #Create list of tuples which contain genotypes of each individual
    GT = zip(*reader)
    f.close()
    sys.stdout.write("Files read\t"+str(time.clock())+"\n")
    #Remove first tuple of position names
    del GT[0]

    #Haplotype length
    hlength=len(GT[0])-1

#Initialize empty tree
G = Tree(0)

#Randomly impute and phase data for first input into tree
for i in GT:
    gt = []
    for n, j in enumerate(i):
        if j[1] == "|":
            a = j[0]
            b = j[2]                
        elif j[1] == "/":
            rand = random.randrange(0,2)
            if rand == 0:
                a=j[0]
                b=j[2]
            else:
                a=j[2]
                b=j[0]         
        else:
            print j[1]
            print "Unknown symbol"

        if a == '.':
            a = alleles[n][random_weighted_choice(allelefreq[n])]
        if b == '.':
            b = alleles[n][random_weighted_choice(allelefreq[n])]
        
        gt.append([a,b])
          
    G.add_genotype(gt)

G = G.merge()

iterations = 10

sys.stdout.write("1st Tree Built\t"+str(time.clock())+"\n")


def sequence(G):
    if G.d == 0:        
        H = Tree(1)
        for i in GT:
            H.add_genotype(forwardbackward(G,i))
    elif G.d == 1:
        H = Tree(0)
        for i in GT:
            j = i[::-1]
            H.add_genotype(forwardbackward(G,j))
        
    else:
        print "Error in graph direction"
    H.merge()
    
    return H
    
while iterations > 0:
    G = sequence(G) 
    sys.stdout.write(str(iterations)+"\t"+str(time.clock())+"\n")

        
    iterations -= 1

output = []
if G.d == 0:
    for i in GT:
        output.append(viterbi(G,i))
    
if G.d == 1:
    for i in GT:        
        j = i[::-1]
        output.append(viterbi(G,j)[::-1])       



for i in output:
    for j, k in enumerate(i):        
        imputed[j] += (str(k[0])+"|"+str(k[1])+"\t")
        

    
g = open(args.output+".txt","w")   


for i in imputed:
    g.write(i+"\n")
     
g.close()
    
sys.stdout.write("End\t"+str(time.clock())+"\n")


