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
        print l
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
        

    #Create list of nodes on each level.
    gnodes = [[]]*(hlength+2)
    

    #Create networkx MultiGraph
    G=nx.MultiDiGraph()

    #Add start node 1 and split the first node.
    G.add_node(1,hap=h,level=0)
    gnodes[0].append(1)
    print gnodes
    print gnodes.append(2)
    print gnodes
    nodesplit(G,1,gnodes)

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
