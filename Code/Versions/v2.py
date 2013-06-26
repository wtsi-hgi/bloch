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

#Function which creates simple sub-tree from node to be used in mergetest
def createtree(n):
    T=nx.DiGraph()
    slevel=[1]*(len(G.node[n]['haplotype'][0])+1) 
    T.add_node(1)    

    for i, hap in enumerate(G.node[n]['haplotype']):
        cnode = 1
       
        for j, h in enumerate(hap):
            
            if h in [(edata['allele']) for u,v,edata in T.out_edges(cnode,data=True) if 'allele' in edata]:
                T[cnode][v]['weight'] = T[cnode][v]['weight'] + G.node[n]['frequency'][i]
                
                cnode = v
            else:
                     
                if j == len(hap)-1 or slevel[j]==slevel[j+1]:
                    pass
                else:
                    #Define function to relabel node names
                    def mapping(x):
                        if x > slevel[j+1]:
                            return x+1
                        else:
                            return x
            
                    #Relabel node names which are above the level of the node being added.
                    nx.relabel_nodes(T,mapping,copy=False)            

                T.add_edge(cnode, slevel[j+1]+1, weight=G.node[n]['frequency'][i], allele=h)
                cnode = slevel[j+1]+1
                for a in range(j+1,len(slevel)):
                    slevel[a] = slevel[a] + 1
    return T
    

#Function which tests if two nodes are similar enough to merge. Returns True or False.
def mergetest(a,b):
    print "mergetest carried out for nodes:"
    print a,b

    #Create 2 trees one, beginning at node a and one beginning at node b.
    A = createtree(a)
    B = createtree(b)

    print "Tree A"
    for i in  A.edges(data=True):
        print i
    print "Tree B"
    for i in B.edges(data=True):
        print i
   
          
    #Score threshold and haplotype count for nodes a and b are calculated.
    totala = sum(G.node[a]['frequency'])
    totalb = sum(G.node[b]['frequency'])

    
    print "totala, totalb"
    print totala, totalb  
    
   
    if totala == 0 and totalb == 0:
        return False
    
    threshold = math.sqrt(math.pow(totala, -1) + math.pow(totalb, -1))
    #maxscore = 0

    #Queue list is created for breadth-first search.
    #q = [[a,b]]

    #Similarity score is calculated for the first item of the queue as long as there are elements in the queue.
    #while len(q) != 0:

    #score = similarityscore(q[0],totala,totalb,threshold)
        
        #If a pair of nodes are not similar enough mergetest returns False immediately.
        #if score is False:
        #return False

        #else:
        #if score > maxscore:
        #maxscore = score
            
            #Each allele of outgoing edges of first node is compared to alleles from outgoing edges of second node.
            #for i in M.out_edges(q[0][0], data=True):
            #x = i[2]['allele']
            #y = None

            #for j in M.out_edges(q[0][1], data=True):
                    
                    #If matching allele is found in outgoing edges of second node then for loop stops.
                    #if j[2]['allele'] == x:
                    #y = j
                    #break
                                                
                    #else:
                    #continue
                    
                #Corresponding pairs of nodes iwth matching alleles are added to the queue.
                #if y is not None:
                #q.append([i[1],y[1]])

            #Completed pair of nodes is removed from list          
            #q.pop(0)
            
    #If all possible pairs are compared and are True, function returns True.
    #return maxscore

            
def merge(l):
    nn = level[l] - level[l-1]
    print "nn"
    print nn

    if nn == 1:
        print "one node in level so no merge can be performed"
        pass
        

    #If there is only 2 nodes in level. Test similarity between the pair of nodes.
    elif nn == 2:
        print "2 nodes in level"
               
        if mergetest(level[l],level[l]-1) == False:
            print "mergetest is false"
        
        else:
            print "nodes need to be merged"
            #mergenodes(level[i],level[i]-1)
            

        for i in range(level[l-1]+1,level[l]+1):
            print i
    else:
        print "more than two nodes in level"

#Splitting function takes a node as a variable and splits node according to haplotypes attached to the node.
def nodesplit(n,l):
    
    #For each haplotype   
    for i in range(haplonum):
        #If edge exists, append frequency and haplotype suffix to the lists on next node.
        if G.node[n]['haplotype'][i][0] in [(edata['allele']) for u,v,edata in G.out_edges(n,data=True) if 'allele' in edata]:
            G.node[v]['haplotype'].append(G.node[n]['haplotype'][i][1:])
            G.node[v]['frequency'].append(G.node[n]['frequency'][i])
        #If edge does not exist, add edge.
        else:
            G.add_edge(n,level[l+1]+1,allele=G.node[n]['haplotype'][i][0])
            G.add_node(level[l+1]+1,haplotype=[G.node[n]['haplotype'][i][1:]],frequency=[G.node[n]['frequency'][i]])
            for a in range(m,(haplolength+1)):
                level[a] = level[a] + 1

#Split all nodes on given level.
def split(l):
    #If first level is being split, then there is only one node on the first level.
    if l == 0:
        nodesplit(1,l)
    #If any other level being split, split each node on the given level
    else:
        for i in range(level[l-1]+1,level[l]+1):
            nodesplit(i,l)
        
        


haplotype = ['1111', '1112', '1122', '1221', '2111', '2112', '2122']
frequency = [21, 79, 95, 116, 25, 112, 152]
haplolength=4
haplonum=7

level=[1]*(haplolength+1)
G=nx.MultiDiGraph()

#Add start node 1
G.add_node(1,haplotype=haplotype,frequency=frequency)


#Split all nodes in level 0
print "split(0)"
split(0)


print 'Edges in Graph'
for i in G.edges(data=True):
    print i
print 'Nodes in Graph'
for i in G.nodes_iter(data=True):
    print i
print 'level'
print level



#Carry out merge part on all nodes in level 1
print "merge(1)"
merge(1)



