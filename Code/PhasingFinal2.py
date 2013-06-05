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

def treealgorithm(h):

    #Function which takes node and splits node according to haplotype dictionary attached to the node
    def nodesplit(G,n,nodes):
        l = G.node[n]['level']
        
        #Loop over keys in node's haplotype dictionary
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

        #Remove haplotype information of node that has been split
        del G.node[n]['hap']


    #Functinon tests if two nodes are similar enough to merge. Returns similarity score or false.
    def mergetest(a,b,l):
        #print a, b
        #print G.node[a], G.node[b]
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
            #print q

            #Splits nodes subgraphs
            nodesplit(A,q[0][0],q[0][2],al)
            nodesplit(B,q[0][1],q[0][2],bl)

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
                    b = [v for u,v,edata in B.out_edges(q[0][1],data=True) if edata['a']==i]
                   
                    
                    #Adds next nodes to list if edges has passed edge test            
                    if len(random.choice(A.node[a[0]]['hap'].keys())) != 0:
                        
                        q.append([a[0],b[0],q[0][2]+1])

            q.pop(0)
          
        return maxs


    #Function merges 2 nodes. a should always be < b
    def mergenodes(a,b):
        print a, b
  
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
                    print "nodes are merged"
                    print minj, mink
                    
                    mergenodes(minj,mink)
                    for a in G.nodes(data=True):
                        print a

                    for a in G.edges(data=True, keys=True):
                        print a
                    #Whenever a merge has taken place, merge is set to True so that the loop is repeated
                    merge = True


    
        
    gnodes = []
    #Create list of nodes on each level.
    for i in range(hlength+2):
        gnodes.append([])
    print gnodes

    #Create networkx MultiGraph
    G=nx.MultiDiGraph()

    #Add start node 1 and split the first node.
    G.add_node(1,hap=h,level=0)
    
    gnodes[0].append(1)
    
    
    nodesplit(G,1,gnodes)
    merge(1)

    
    for i in G.nodes(data=True):
        print i
    for i in G.edges(data=True, keys=True):
        print i
    print gnodes

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
print haplotypes
print treealgorithm(haplotypes)
