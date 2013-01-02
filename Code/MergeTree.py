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
from networkx.utils import is_string_like

from CreateTree import *

#Function which returns the level of the given node
def whichlevel(n):    
    for i in range(len(mlevel)):
       if mlevel[i] <= n <= mlevel[i+1]:
	    return i+1
        
#Function which removes node from nlevel
def remove_node_from_nlevel(m):    
    for a in range(m,(haplolength++1)):
        nlevel[a] = nlevel[a]-1

#Inner merge function
def merge_inner(a,b):
	
    #List of edges directed out of nodes a and b.
    a_edges = M.out_edges(a, data=True)
    b_edges = M.out_edges(b, data=True)

    #Iterate through outgoing edges of node a
    for i in a_edges:
        x = i[2]['allele']
        y = None
	#If outgoing edge of a has same allele as outgoing edge of b the index is set to y.
        for iterator in b_edges:
            if iterator[2]['allele'] == x:
                y = iterator

	#If y exists, it is removed from list of outgoing edges from b
        if y is not None:
            b_edges.remove(y)         
                       
            #get edge from b corresponding to i
            j = y
	    
            #Weight of edge in subtree of b is added to corresponding edge of a
            M.edge[i[0]][i[1]][0]['weight'] = M.edge[i[0]][i[1]][0]['weight'] + M.edge[j[0]][j[1]][0]['weight']
            
            #Recurse using child node
            merge_inner(i[1],j[1])

            #Remove edge and leaf node that is added to subgraph of a.
            M.remove_edge(y[0],y[1])
            M.remove_node(y[1])
            l=whichlevel(y[1])
            remove_node_from_nlevel(l) 

    #All outgoing edges of b which match one from a have been removing leaving only missing edges in list b_edges
    for k in b_edges:

	#Missing edges added to subgraph of a and remove the missing edge from graph
        M.add_weighted_edges_from([(a, k[1], k[2]['weight'])], allele = k[2]['allele'])
        M.remove_edge(k[0],k[1])
   
    
def mergenodes(a, b):
    #Call merge_inner
    merge_inner(a,b)

    #List of edges going into a
    b_in = M.in_edges(b, data=True)

    #Move all incoming edges of b in a.
    for i in b_in:
        M.add_weighted_edges_from([(i[0], a, i[2]['weight'])], allele = i[2]['allele'])
        M.remove_edge(i[0],i[1])

    #Remove node b	
    M.remove_node(b)
    l=whichlevel(b)
    remove_node_from_nlevel(l)
               
    #Relabel node names which are above the level of the node being added.
    mapping=dict(zip(M.nodes(),range(1,nlevel[haplolength]+1)))
    nx.relabel_nodes(M,mapping,copy=False)

      
#Set mlevel to equal level
mlevel = list(level)
nlevel = list(mlevel)


#Create copy of G into multidirected graph
M = nx.MultiDiGraph()
M.add_nodes_from(G)
M.add_edges_from(G.edges_iter(data=True))

 
#Iterate through levels of M starting at second level.
for i in range(1, haplolength+1):

    nn = mlevel[i] - mlevel[i-1]
    
    #If there is only 1 node in graph G level i then carry on to next level. 
    if nn == 1:
      continue

    #If there is only 2 nodes in level. Test similarity between the pair of nodes.
    elif nn == 2:

        if mergetest(mlevel[i],mlevel[i]-1) == False:
            continue
        
        else:
            mergenodes(mlevel[i],mlevel[i-1])       
            

    #If there is more than one node in level, test each pair of nodes in each level and merge 
    else:


