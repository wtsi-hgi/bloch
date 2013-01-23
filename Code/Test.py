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
    a_out = M.out_edges(a, data=True)
    b_out = M.out_edges(b, data=True)

    #Iterate through outgoing edges of node a
    for i in a_out:
        x = i[2]['allele']
        y = None
        
	#If outgoing edge of a has same allele as outgoing edge of b the index is set to y.
        for iterator in b_out:
            if iterator[2]['allele'] == x:
                y = iterator

	#If y exists, it is removed from list of outgoing edges from b
        if y is not None:
            b_out.remove(y)         
                       
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

    #All outgoing edges of b which match one from a have been removing leaving only missing edges in list b_out
    for k in b_out:

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
    print 'reached this point'
    #Remove node b	
    M.remove_node(b)
    l=whichlevel(b)
    remove_node_from_nlevel(l)
               
    #Relabel node names which are above the level of the node being added.
    mapping=dict(zip(M.nodes(),range(1,nlevel[haplolength]+1)))
    nx.relabel_nodes(M,mapping,copy=False)
    
    
    
  
def plot_tree(G, name='plot.png'):
      
    #Create undirected graph from directed graph for visualization purposes           
    N=G.to_undirected()
    pyl.rcParams['figure.figsize'] = 10, 10
    
    #Clear plot to begin again
    plt.clf()
    plt.axis('off')
        
    #Set vizualization layout
    prog ='dot'
    pos = nx.drawing.graphviz_layout(G, prog)
    
    #Mark edges with labels corresponding to the weight of the edge
    edge_labels=dict([((u,v,),d['weight']) for u,v,d in N.edges(data=True)])
   
    #Draw out nodes
    nx.draw_networkx_nodes(G,pos,node_size=400, node_color='w')
    #Draw edge labels
    nx.draw_networkx_edge_labels(G,pos,ax=None,edge_labels=edge_labels)
    
    #Create two edges lists for the different alleles
    allele_1 = [(u,v) for (u,v,d) in N.edges(data=True) if d['allele'] == 1]
    allele_2 = [(u,v) for (u,v,d) in N.edges(data=True) if d['allele'] == 2]
    #both_alleles = [(u,v) for (u,v,d) in N.edges(data=True) if d['weight'] >0.5]
    #Draw node labels
    nx.draw_networkx_labels(G,pos,ax=None)

    #Draw edges from edge lists of different alleles
    nx.draw_networkx_edges(N, pos, edgelist=allele_1, node_size=100, width=3, with_labels=False,ax=None)
    nx.draw_networkx_edges(N, pos, edgelist=allele_2, node_size=100, width=6, with_labels=False, style='dashed',ax=None)
    #nx.draw(N, pos, node_size=100, node_color='w', edge_color=edge_colours, width=4, with_labels=False)

    #os.remove('plot.png')
    #Show plot in window
    #plt.show()

    #Save plot
    plt.savefig(name)
    
	    
#Create test graph manually
M = nx.MultiDiGraph()
M.add_weighted_edges_from([(1,2,3.5),(2,4,5.0),(4,8,9.0),(5,10,2.7),(3,6,6.0),(6,11,7.0),(7,12,4.0)], allele=1)
M.add_weighted_edges_from([(1,3,5.0),(2,5,1.0),(4,9,2.0),(3,7,8.0),(7,13,1.0)], allele=2)

#Set mlevel to
mlevel = [1,3,7,13]
nlevel=list(mlevel)
haplolength = 3

plot_tree(M, name='before_merge.png')
#plot_pgv(A)

mergenodes(2,3)

plot_tree(M, name='after_merge.png')

O=nx.to_agraph(A)
O.draw('plot.png', prog='nop')

plot_tree(A)

