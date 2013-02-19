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
import networkx as nx
import sys
import matplotlib.pyplot as plt
import pygraphviz as pgv
import operator
import os
import math

#from CheckInput import *

#Splitting function takes a node as a variable and splits node according to haplotypes attached to the node.
def nodesplit(G,n,l,level):
    
    #For each haplotype   
    for i in range(len(G.node[n]['haplotype'])):
        #If edge exists, append frequency and haplotype suffix to the lists on next node.
        if G.node[n]['haplotype'][i][0] in [(edata['allele']) for u,v,edata in G.out_edges(n,data=True) if 'allele' in edata]:
            G.node[v]['haplotype'].append(G.node[n]['haplotype'][i][1:])
            G.node[v]['frequency'].append(G.node[n]['frequency'][i])

        #If edge does not exist, add edge.
        else:
            G.add_edge(n,level[l+1]+1,allele=G.node[n]['haplotype'][i][0])
            G.add_node(level[l+1]+1,haplotype=[G.node[n]['haplotype'][i][1:]],frequency=[G.node[n]['frequency'][i]])

            for a in range(l+1,(len(G.node[n]['haplotype'][0])+1)+l):
                level[a] = level[a] + 1      

#Function which tests if two nodes are similar enough to merge. Returns True or False.
def mergetest(a,b):
    #Score threshold and haplotype count for nodes a and b are calculated.
    totala = sum(G.node[a]['frequency'])
    totalb = sum(G.node[b]['frequency'])
    threshold = math.sqrt(math.pow(totala, -1) + math.pow(totalb, -1))
    maxscore = 0
    alevel = [1]*(len(G.node[a]['haplotype'][0])+1)
    blevel = [1]*(len(G.node[a]['haplotype'][0])+1)
    
    A=nx.DiGraph()
    A.clear()
    A.add_node(1,haplotype=G.node[a]['haplotype'],frequency=G.node[a]['frequency'])
    B=nx.DiGraph()
    B.clear()
    B.add_node(1,haplotype=G.node[b]['haplotype'],frequency=G.node[b]['frequency'])
        
    q = [[1,1,0]]

    while len(q) != 0:
        nodesplit(A,q[0][0],q[0][2],alevel)
        nodesplit(B,q[0][1],q[0][2],blevel)
        
        dicta = dict((edata['allele'],sum(A.node[v]['frequency'])) for u,v,edata in A.out_edges(q[0][0], data=True))
        dictb = dict((edata['allele'],sum(B.node[v]['frequency'])) for u,v,edata in B.out_edges(q[0][0], data=True))
        
        seta = set(dicta.keys())
        setb = set(dictb.keys())
        
        for i in list(seta.symmetric_difference(setb)):
            if i in dicta:

                score = math.fabs(float(dicta[i])/float(totala))
                if score > threshold:

                    return False
                elif score > maxscore:
                    maxscore = score

            if i in dictb:

                score = math.fabs(float(dictb[i])/float(totalb))
                if score > threshold:
    
                    return False
                elif score > maxscore:
                    maxscore = score
                   
        for i in list(seta.intersection(setb)):
            score = math.fabs(float(dicta[i])/totala - (float(dictb[i])/totalb))
            if score > threshold:   
                return False
                
            else:
                if score > maxscore:
                    maxscore = score          
               
                a = [v for u,v,edata in A.out_edges(q[0][0],data=True) if edata['allele']==i]
                               
                b = [v for u,v,edata in A.out_edges(q[0][0],data=True) if edata['allele']==i]
                                  
                if len(A.node[a[0]]['haplotype'][0]) != 0:
                    q.append([a[0] ,b[0] ,q[0][2]+1])

        q.pop(0)
  
    return maxscore

def mergenodes(a,b):
    for i in range(len(G.node[b]['haplotype'])):
          
        if G.node[b]['haplotype'][i] in [j for j in G.node[a]['haplotype']]:
            index = [j for j in G.node[a]['haplotype']].index(G.node[b]['haplotype'][i])
            G.node[a]['frequency'][index] = G.node[a]['frequency'][index] + G.node[b]['frequency'][i]

        else:
            G.node[a]['haplotype'].append(G.node[b]['haplotype'][i])
            G.node[a]['frequency'].append(G.node[b]['frequency'][i])

    for i in G.in_edges(b, data=True):
        G.add_edge(i[0],a,allele=i[2]['allele'])
        G.remove_edge(i[0],i[1])

    G.remove_node(b)

    #Edit glevel!!!
    for i in glevel:
        if b <= i:
            glevel[glevel.index(i)] = glevel[glevel.index(i)] - 1

def merge(l):
    nn = glevel[l] - glevel[l-1]
    if nn == 1:     
        pass
    
    #If there is only 2 nodes in level. Test similarity between the pair of nodes.
    elif nn == 2:
               
        if mergetest(glevel[l],glevel[l]-1) == False:
            pass
        
        else:
            mergenodes(level[i],level[i]-1)

    else:
        merge = True
        while merge == True:
            merge = False
            levelmin = 1000
            minj = 0
            mink = 0

            #Calculate merge score for each pair of nodes
            for j in range(glevel[l-1]+1, glevel[l]+1):
                for k in range(j, glevel[l]+1):               
                    
                    if j != k:
                        testscore = mergetest(j,k)

                 
                        if testscore == False:
                            continue
                        else:
                            if testscore < levelmin:
                                levelmin = testscore
                                minj = j
                                mink = k
                                            
            if levelmin != 1000:
                mergenodes(minj,mink)
                merge = True

haplotype = ['1111', '1112', '1122', '1221', '2111', '2112', '2122']
frequency = [21, 79, 95, 116, 25, 112, 152]
haplolength=4
haplonum=7

glevel=[1]*(haplolength+1)
G=nx.MultiDiGraph()

#Add start node 1
G.add_node(1,haplotype=haplotype,frequency=frequency)
nodesplit(G,1,0,glevel)

for i in range(1,haplolength):
    merge(i)

    for j in range(glevel[i-1]+1,glevel[i]+1):
        nodesplit(G,j,i,glevel)
    
print 'Edges in Graph'
for i in G.edges(data=True):
    print i
print 'Nodes in Graph'
for i in G.nodes_iter(data=True):
    print i
print 'glevel'
print glevel
