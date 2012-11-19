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

from CheckInput import *

#Create list of levels. Each position is a cumulative count of the number
#of nodes on the level and the previous levels.
level = [0]*(haplolength+1)

#Set current node pointer
cnode = 0

#Create empty directed graph
G = nx.DiGraph()

#Function to add node to level count. l=level
def add_node_to_level(m):
    
    for a in range(m,(haplolength+1)):
        level[a] = level[a] + 1      
        
#Function to add node to current node c=cnode, h=haplotype, l=level
def add_edge_from_cnode(G,c,h,l):
    
    #If node is being created at last level
    if l == haplolength - 1:
        
        #Add edge from current node to node with name one higher than existing nodes.
        #Set weight and allele to edge data.
        G.add_edge(c, (level[l+1]+1), weight=frequency[h], allele=haplotype[h][l])
        
    #If node is being created that is not at the last level
    else:
        
        #Define function to relabel node names
        def mapping(x):
            if x > level[l+1]:
                return x+1
            else:
                return x
            
        #Relabel node names which are above the level of the node being added.
        nx.relabel_nodes(G,mapping,copy=False)
            
        #Add edge from current node to node with name one higher than the last node on the next level.
        #Set weight and allele to edge data.
        G.add_edge(c, (level[l+1]+1), weight=frequency[h], allele=haplotype[h][l])

#Function which checks if edge from node exists. n=node to be checked. a=allele.
def check_for_existing_edge(n, a):

    #Iterates through successor nodes of node n
    for i in G.successors_iter(n):
        #If allele matches successor node then the sucessor node name is returned
        if G[n][i]['allele'] == a:            
            return i
    #If no successor nodes match the allele then function returns 0 
    return 0


#For each haplotype
for i in range(haplonum):
    cnode = 1

    #For each allele in haplotype
    for j in range(haplolength):
                       
        #If no nodes have been created, add first edge between node 1 and node 2
        if level[haplolength] == 0:
            G.add_edge(1, 2, weight=frequency[i], allele=haplotype[i][j])
            add_node_to_level(0)
            add_node_to_level(1)
            cnode = 2

        #If current node has no outgoing edges and an edge between current node and next level.
        elif len(G.successors(cnode)) == 0:            
            add_edge_from_cnode(G,cnode,i, j)
            cnode = level[j+1]+1
            add_node_to_level(j+1)      
                                
        #If current node has at least one outgoing edge 
        else:
            
            #Check if there is an existing edge with the same allele value.
            snode = check_for_existing_edge(cnode, haplotype[i][j])                 
                
            #If no edges have the same allele and an edge between current node and next level.
            if snode == 0:           
                add_edge_from_cnode(G,cnode,i, j)
                add_node_to_level(j+1)
                cnode = level[j+1]
                    
            #If there exists an edge with the same allele add corresponding weight to edge.
            else:
                G[cnode][snode]['weight'] = G[cnode][snode]['weight'] + frequency[i]
                cnode = snode

#Create undirected graph from directed graph for visualization purposes           
H=G.to_undirected()        
                    
#Set vizualization layout
prog='dot'
pos=nx.drawing.graphviz_layout(G, prog)

#Mark edges with labels corresponding to the weight of the edge
edge_labels=dict([((u,v,),d['weight']) for u,v,d in H.edges(data=True)])

#Draw edge labelsquit()

nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_labels)

#Create edge_colours list and iterate through edges of graph H
edge_colours=[]
for u,v,d in H.edges(data=True):
    #If allele is labelled 1 append colour MidnightBlue
    if d['allele'] == '1':
        edge_colours.append('MidnightBlue')
        #Otherwise append the colour LightBlue
    else:
        edge_colours.append('LightBlue')

#Draw the graph using edge_colours calculated previously
nx.draw(H, pos, node_size=100, node_color='w', edge_color=edge_colours, width=4, with_labels=False)


#os.remove('plot.png')

#plt.show()
#Show plot in window
#plt.savefig("plot.png")

print "level" + str(level)

print "dfs_edges"
for i in nx.dfs_edges(G):
    print i




print "nx.bfs_edges(G, 1)"
for i in nx.bfs_edges(G, 1):
    print i

hapnum = len(list(set(haplotype)))
print hapnum

