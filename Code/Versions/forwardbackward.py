#Import required python packages
import csv
import networkx as nx
import random
import numpy as np
import itertools
import math

from PhasingFinal import T

gt = ('.|.', '0|0', '0|1', '0|1')

#Haploid initial state probabilities. allele: allele number.
def hapinitial(allele):    
    #Edge counts are used. Count of edge/Total count for all edges.
    for i in T[0].out_edges(1, data=True):        
        if i[2]['a'] == allele:           
            return float(i[2]['weight'])/float(sum(T[0].node[1]['hap'].itervalues()))

#Diploid initial state probabilities. a,b: allele number
def dipinitial(a,b):    
    return hapinitial(a)*hapinitial(b)

#Emission state probabilities. gt: genotype. s: tuple of alleles.
def emission(gt,s):
    #If one allele is unknown. If known allele is contained in s then 1 is returned.
    if gt.count('.') == 1:
        for i in gt:
            if (i == "|") or (i == "\\"):
                continue
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
        if set([gt[0],gt[-1]]) == set(s):
            return 1
        else:
            return 0

#Transition state probabilities. e,d: edge tuples incl data
def haptrans((e,d)):
    #If parent node of edge e is child node of edge d.
    if e[0] == d[1]:
            
        #Edge count/parent node count is returned
        return float(e[3]['weight'])/float(sum(T[0].node[e[0]]['hap'].itervalues())) 
    else:
        return 0.0

#Diploid transition probabilities
def diptrans(a,b):
    return haptrans(a)*haptrans(b)

#Function to find edge corresponding to matrix coordinate. g=flattened index. l=level(index of m)
def findedge(g,l):    
    #Convert flattened index to coordinate
    t = np.unravel_index(g, (n[l],n[l]))    
    f = [-1,-1]
    if l == 0:
        x = 1
    else:
        x = [j for j in range(T[1][l-1]+1,T[1][l]+1)]        
        
    #Set f to equal tuple of ordered edges that the index g corresponds to
    for a, b in enumerate(T[0].out_edges(nbunch=x, keys=True, data=True)):
        if t[0] == a:
                f[0] = b
        if t[1] == a:
                f[1] = b

    if -1 in f:
        print 'Error in findedge'
    else:
        return f

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
    if emission(gt[0],(a[1][3]['a'],b[1][3]['a'])) != 0:              
        if a[1] == b[1]:
            t = hapinitial(a[1][3]['a'])
            var = t*t            
        else:
           var = dipinitial(a[1][3]['a'], b[1][3]['a'])

        #Matrix element is set to calculated diploid initial probability
        m[0][a[0]][b[0]] = var
        

#Induction. Iterate through each level
for i in range(1,hlength+1):

    #Iterate through ordered pairs of outgoing edges in each level
    for a, b in itertools.product([(c, d) for c, d in enumerate(T[0].out_edges(nbunch=[j for j in range(T[1][i-1]+1,T[1][i]+1)], keys=True, data=True))], repeat=2):
        
        #If emission probability does not equal 0
        if emission(gt[i],(a[1][3]['a'],b[1][3]['a'])) != 0.0:
            print "passed emission test"
            
            if i == 1:
                nodes = 1
            else:                
                nodes = [j for j in range(T[1][i-2]+1,T[1][i-1]+1)]
                
            var = sum([(m[i-1][c[0]][d[0]]*diptrans([a[1],c[1]], [b[1],d[1]])) for c, d  in itertools.product([(e, f) for e, f in enumerate(T[0].out_edges(nbunch=nodes, keys=True, data=True))], repeat=2)])
            #Matrix element is set to var calculation formula
            m[i][a[0]][b[0]] = var

for i in m:
    print i
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
    
    e = findedge(s[i], i)

    #For each position in forward probability matrix of level below, calculate sampling probabilities
    for b, j in enumerate(m[i-1].flat):
        d = findedge(b, i-1)
        prob.append((emission(gt[i],(e[0][3]['a'],e[1][3]['a']))*diptrans((e[0],d[0]),(e[1],d[1]))*j)/m[i].flat[s[i]])

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
    e = findedge(s[i],i)
    sample[0] += e[0][3]['a']
    sample[1] += e[1][3]['a']

return sample
