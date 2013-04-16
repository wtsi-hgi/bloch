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
from __future__ import division
import networkx as nx
import sys
import operator
import os
import math
import numpy as np
import itertools

from v3 import G, glevel

#ll is number of levels in the graph
ll = len(glevel) -1

#Set the genotype
GT = (('?','?'),(1,1),(1,2),(1,2))

#Create n the list of the number of edges in each level
n = [G.out_degree(1)]

#m is the list of matrices of forward probabilities at each level
m = [np.zeros(shape=(n1,n1))]

#Append matrices to list m
for i in range(ll-1):
    
    #Sum n for each level
    n.append(G.out_degree(glevel[i]+1))
    for j in range(glevel[i]+2, glevel[i+1]+1):
        n[i+1] += G.out_degree(j)

    #Append zero-filled matrix to list m
    m.append(np.zeros(shape=(n[i+1],n[i+1])))

#Haploid initial state probabilities
def hapinitial(allele):
    
    #Edge counts are used. Count of edge/Total count for all edges.
    for i in G.out_edges(1, data=True):        
        if i[2]['allele'] == allele:
            return sum(G.node[i[1]]['frequency'])/sum(G.node[1]['frequency'])

#Diploid initial state probabilities
def dipinitial(a,b):    
    return hapinitial(a)*hapinitial(b)

#Emission state probabilities
def emission(gt,s):
    if gt.count('?') == 1:
        for i in gt:
            if i != '?':
                if i in s:
                    return 1
                else:
                    return 0        
    elif gt.count('?') == 2:        
        return 1    
    else:
        if set(gt) == set(s):
            return 1
        else:
            return 0

#Transition state probabilities
def haptrans((e,d)):
    if e[0] == d[1]: #parent node of edge e is child node of edge d
        if 'weight' in e[3]:
            x = e[3]['weight']
        else:
            x = sum(G.node[e[1]]['frequency'])            
        return x/sum(G.node[e[0]]['frequency']) 
    else:
        return 0.0

#Diploid transition probabilities
def diptrans(a,b):
    return haptrans(a)*haptrans(b)

#Initiation
for i in G[1].keys():
    for j in G[1][i].keys():
        print G[1][i][j]

for a, b  in itertools.product([(i, j) for i, j in enumerate(G.out_edges(1, keys=True, data=True))], repeat=2):
    print (a[1],b[1])
    if emission(GT[0],(a[1][3]['allele'],b[1][3]['allele'])) == 1:        
        if a[1] == b[1]:
            t = hapinitial(a[1][3]['allele'])
            var = t*t            
        else:
            var = dipinitial(a[1][3]['allele'], b[1][3]['allele'])
        a1[a[0]][b[0]]=var

#Induction
#Iterate through each level
for i in range(1,ll):

#Iterate through ordered pairs of outgoing edges in each level
    for a, b in itertools.product([(c, d) for c, d in enumerate(G.out_edges(nbunch=[j for j in range(glevel[i-1]+1,glevel[i]+1)], keys=True, data=True))], repeat=2):
        
        #Check if emission probability = 1
        if emission(GT[i],(a[1][3]['allele'],b[1][3]['allele'])) == 1:
            
            if i == 1:
                nodes = 1
            else:
                nodes = [j for j in range(glevel[i-2]+1,glevel[i-1]+1)]

            var = sum([(globals()['a'+str(i)][c[0]][d[0]]*diptrans((a[1],c[1]), (b[1],d[1]))) for c, d  in itertools.product([(e, f) for e, f in enumerate(G.out_edges(nbunch=nodes, keys=True, data=True))], repeat=2)])

            globals()['a'+str(i+1)][a[0]][b[0]] = var


print m[0]
print m[1]
print m[2]
print m[3]


