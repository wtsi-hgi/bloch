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
        
        

#Inner merge function
def merge_inner(a,b):
	
    #List of edges directed out of nodes a and b.
    a_edges = A.out_edges(a, data=True)
    b_edges = A.out_edges(b, data=True)

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
            A.edge[i[0]][i[1]][0]['weight'] = A.edge[i[0]][i[1]][0]['weight'] + A.edge[j[0]][j[1]][0]['weight']
            
            #Recurse using child node
            merge_inner(i[1],j[1])

           
	            

       

    #All outgoing edges of b which match one from a have been removing leaving only missing edges in list b_edges
    for k in b_edges:
	    
	#Missing edges added to subgraph of a and remove the missing edge from graph
        A.add_weighted_edges_from([(a, k[1], k[2]['weight'])], allele = k[2]['allele'])
        A.remove_edge(k[0],k[1])
   
    
def merge(a, b):
    #Call merge_inner
    merge_inner(a,b)

    #List of edges going into a
    b_in = A.in_edges(b, data=True)

    #Move all incoming edges of b in a.
    for i in b_in:
        A.add_weighted_edges_from([(i[0], a, i[2]['weight'])], allele = i[2]['allele'])
        A.remove_edge(i[0],i[1])

    #Remove node b	
    A.remove_node(b)
    l=whichlevel(b)
            
    remove_node_from_nlevel(l)


    #Rename node lables so that they are in number order and have no gaps. edit mlevel list.
    print mlevel
    print nlevel
    print A.nodes()

   
    mapping=dict(zip(A.nodes(),range(1,nlevel[haplolength]+1)))
            
    #Relabel node names which are above the level of the node being added.
    nx.relabel_nodes(A,mapping,copy=False)
    print A.nodes()
    
    
    
  
def plot_tree(G, name='plot.png'):
    
    
    #Create undirected graph from directed graph for visualization purposes           
    N=G.to_undirected()        

    #Clear plot to begin again
    plt.clf()
    
    #Set vizualization layout
    prog ='dot'
    pos = nx.drawing.graphviz_layout(G, prog)

    #Mark edges with labels corresponding to the weight of the edge
    edge_labels=dict([((u,v,),d['weight']) for u,v,d in N.edges(data=True)])

    #Draw edge labels
    nx.draw_networkx_edge_labels(G,pos,ax=None,edge_labels=edge_labels)

    #Create two edges lists for the different alleles
    allele_1 = [(u,v) for (u,v,d) in N.edges(data=True) if d['allele'] == 1]
    allele_2 = [(u,v) for (u,v,d) in N.edges(data=True) if d['allele'] == 2]
    #both_alleles = [(u,v) for (u,v,d) in N.edges(data=True) if d['weight'] >0.5]
    
    #Draw out nodes
    nx.draw_networkx_nodes(G,pos,node_size=500, node_color='w',ax=None)

    #Draw node labels
    nx.draw_networkx_labels(G,pos,ax=None)

    #Draw edges from edge lists of different alleles
    nx.draw_networkx_edges(N, pos, edgelist=allele_1, node_size=100, width=2, with_labels=False,ax=None)
    nx.draw_networkx_edges(N, pos, edgelist=allele_2, node_size=100, width=6, with_labels=False, style='dashed',ax=None)
    #nx.draw(N, pos, node_size=100, node_color='w', edge_color=edge_colours, width=4, with_labels=False)


    #os.remove('plot.png')
    #Show plot in window
    #plt.show()

    #Save plot
    plt.savefig(name)

    #def whichlevel(n):
    #for i in mlevel:


    #def plot_pgv(G, name='plot.png'):
    
	    
#Create test graph manually
A = nx.MultiDiGraph()
A.add_weighted_edges_from([(1,2,3.5),(2,4,5.0),(4,8,9.0),(5,10,2.7),(3,6,6.0),(6,11,7.0),(7,12,4.0)], allele=1)
A.add_weighted_edges_from([(1,3,5.0),(2,5,1.0),(4,9,2.0),(3,7,8.0),(7,13,1.0)], allele=2)

#Set mlevel to
mlevel = [1,3,7,13]
nlevel=list(mlevel)
haplolength = 3

plot_tree(A, name='before_merge.png')
#plot_pgv(A)

merge(2,3)

plot_tree(A, name='after_merge.png')

O=nx.to_agraph(A)
O.draw('plot.png', prog='nop')









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




plot_tree(A)

