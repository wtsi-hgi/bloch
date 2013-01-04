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
import pylab as pyl
import pygraphviz as pgv
import operator
import os
import math
from networkx.utils import is_string_like
from collections import deque


def mergetest(a,b,allele):
    totala = math.fsum([i[2]['weight'] for i in A.out_edges(a, data=True)])
    totalb = math.fsum([i[2]['weight'] for i in A.out_edges(b, data=True)])
    threshold = math.sqrt(math.pow(totala, -1) + math.pow(totalb, -1))

    for i in A.out_edges(a, data=True):
        if i[2]['allele'] == allele:
            aw = i[2]['weight']
        else:
            aw = 0

    for i in A.out_edges(b, data=True):
        if i[2]['allele'] == allele:
            bw = i[2]['weight']
        else:
            bw = 0

    score = maths.fabs(aw/totala - bw/totalb)

    if score > threshold:
        return False
    else:
        return True

    
    #Queue root node
    #While queue is not empty
    #Carry out edge test on all decendants
    #If false return false
    #If not false add decendants to queue and delete root nodes etc

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
A = nx.MultiDiGraph()
A.add_weighted_edges_from([(1,2,311),(2,4,195),(4,7,100),(7,12,21),(3,6,289),(6,10,137),(10,16,25)], allele=1)
A.add_weighted_edges_from([(7,13,79),(4,8,95),(8,14,95),(2,5,116),(5,9,116),(9,15,116),(1,3,289),(10,17,112),(6,11,152),(11,18,152)], allele=2)

#Set mlevel to
mlevel = [1,3,6,11,18]
haplolength = 4

plot_tree(A)


