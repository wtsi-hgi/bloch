#Import required python packages
import csv
import networkx as nx
import random
import numpy as np
import itertools
import math

from PhasingFinal import T

gt = ('.|.', '1|1', '1|2', '1|2')

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

#Haploid initial state probabilities. allele: allele number.
def hapinitial(allele):    
    #Edge counts are used. Count of edge/Total count for all edges.
    for i in T[0].out_edges(1, data=True):        
        if i[2]['a'] == allele:
            print T[0].node[i[1]]['hap']
            return sum(T[0].node[i[1]]['hap'].itervalues())/sum(T[0].node[1]['hap'].itervalues())

#Diploid initial state probabilities. a,b: allele number
def dipinitial(a,b):    
    return hapinitial(a)*hapinitial(b)

#Emission state probabilities. gt: genotype. s: tuple of alleles.
def emission(gt,s):    
    #If one allele is unknown. If known allele is contained in s then 1 is returned.
    if gt.count('.') == 1:
        for i in gt:
            if i != '.':
                if i in s:
                    return 1
                else:
                    return 0
                
    #If both alleles are unknown 1 is returned
    elif gt.count('.') == 2:        
        return 1

    #If gt equals s without order 1 is returned
    else:
        if set(gt) == set(s):
            return 1
        else:
            return 0

#Transition state probabilities. e,d: edge tuples incl data
def haptrans((e,d)):
    #If parent node of edge e is child node of edge d.
    if e[0] == d[1]:         
        #Edge count/parent node count is returned
        return e[3]['weight']/sum(T[0].node[e[0]]['hap'].itervalues()) 
    else:
        return 0.0

#Diploid transition probabilities
def diptrans(a,b):
    return haptrans(a)*haptrans(b)

#Initiation. Iterate through pairs of outgoing edges from node 1.
for a, b  in itertools.product([(i, j) for i, j in enumerate(T[0].out_edges(1, keys=True, data=True))], repeat=2):    
    #If emmsion probability does not equal 0
    if emission(gt[0],(a[1][3]['a'],b[1][3]['a'])) != 0:              
        if a[1] == b[1]:
            t = hapinitial(a[1][3]['a'])
            
            var = t*t            
        else:
            var = dipinitial(a[1][3]['a'], b[1][3]['a'])

        #Matrix element is set to calculated diploid initial probability
        m[0][a[0]][b[0]] = var

#Induction. Iterate through each level
for i in range(1,hlength-1):
    
    #Iterate through ordered pairs of outgoing edges in each level
    for a, b in itertools.product([(c, d) for c, d in enumerate(T[0].out_edges(nbunch=[j for j in range(T[1][i-1]+1,T[1][i]+1)], keys=True, data=True))], repeat=2):
        
        #If emission probability does not equal 0
        if emission(gt[i],(a[1][3]['a'],b[1][3]['a'])) != 0:
            
            if i == 1:
                nodes = 1
            else:
                nodes = [j for j in range(T[0][i-2]+1,T[0][i-1]+1)]
                
            var = sum([(m[i-1][c[0]][d[0]]*diptrans([a[1],c[1]], [b[1],d[1]])) for c, d  in itertools.product([(e, f) for e, f in enumerate(G.out_edges(nbunch=nodes, keys=True, data=True))], repeat=2)])
            #Matrix element is set to var calculation formula
            m[i][a[0]][b[0]] = var

for i in m:
    print i



